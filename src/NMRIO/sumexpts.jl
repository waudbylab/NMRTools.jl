#!/usr/bin/env julia

"""
    sumexpts(inputexpts...; out, weights=[])

Sum a collection of Bruker NMR experiments, with optional weighting factors.
Data will be truncated to fit the shortest data file.

# Arguments
- `inputexpts...`: Input experiment names (strings or integers)
- `out`: Output experiment name (string or integer)
- `weights`: Optional vector of weights for each input experiment

# Details
- Experiment names should be either strings with filenames or integers converted to strings
- If weights are supplied, the length must match the number of input experiments
- First input experiment is used as a template for the output experiment
- For 1D experiments, 'fid' files are added; for nD experiments, 'ser' files are added
- Processed data files in pdata/1 are removed, other pdata subdirectories are deleted
- Updates the title file with information about the summed experiments
- Updates the number of scans in the acqus file to be the sum of individual experiments
- Prompts before overwriting any existing experiment
"""
function sumexpts(inputexpts...; out, weights=[])
    # Convert experiment names to strings if needed
    outputexpt_str = typeof(out) <: Integer ? string(out) : out
    inputexpts_str = [typeof(expt) <: Integer ? string(expt) : expt for expt in inputexpts]
    
    println("=== NMR Experiment Summation ===")
    println("Output: $outputexpt_str")
    println("Inputs: $(join(inputexpts_str, ", "))")
    
    # Check if input directories exist
    for expt in inputexpts_str
        if !isdir(expt)
            throw(ArgumentError("Input experiment directory does not exist: $expt"))
        end
    end
    
    # Check if weights are supplied and have correct length
    if !isempty(weights) && length(weights) != length(inputexpts)
        throw(DimensionMismatch("Length of weight vector must match the number of input experiments"))
    end
    
    # Create default weights if not supplied
    if isempty(weights)
        weights = ones(length(inputexpts))
    end
    
    println("Weights: $(join(weights, ", "))")
    
    # Check if the output experiment already exists and prompt before overwriting
    if isdir(outputexpt_str)
        print("Output experiment directory already exists. Overwrite? (y/n): ")
        response = lowercase(readline())
        if response != "y" && response != "yes"
            println("Operation cancelled.")
            return nothing
        end
    end
    
    # Copy the first input experiment to create the output experiment
    if isdir(outputexpt_str)
        # Remove the existing output directory first
        rm(outputexpt_str, recursive=true)
        println("Removed existing output experiment: $outputexpt_str")
    end
    
    # Copy the entire directory structure from the first input experiment
    cp(inputexpts_str[1], outputexpt_str)
    _makewritable(outputexpt_str)
    println("Created output experiment from $(inputexpts_str[1])")
    
    # Read binary format parameters from the first experiment's acqus file
    acqus_path = joinpath(inputexpts_str[1], "acqus")
    bytorda = 0  # default: little endian
    dtypa = 0    # default: int32
    
    if isfile(acqus_path)
        open(acqus_path, "r") do f
            for line in eachline(f)
                if occursin("##\$BYTORDA=", line)
                    bytorda = parse(Int, match(r"##\$BYTORDA=\s*(\d+)", line).captures[1])
                elseif occursin("##\$DTYPA=", line)
                    dtypa = parse(Int, match(r"##\$DTYPA=\s*(\d+)", line).captures[1])
                end
            end
        end
    end
    
    # Determine data type based on DTYPA parameter
    datatype = dtypa == 0 ? Int32 : Float64
    endianness = bytorda == 1 ? "big" : "little"
    
    println("Data format: $(datatype) ($(endianness)-endian)")
    
    # Determine if experiment is 1D or nD by checking for existence of 'fid' file
    is_1d = isfile(joinpath(inputexpts_str[1], "fid"))
    data_file = is_1d ? "fid" : "ser"
    println("Detected $(is_1d ? "1D" : "nD") experiment, using $data_file files")
    
    # Clean up pdata directories
    if isdir(joinpath(outputexpt_str, "pdata"))
        # Keep pdata/1 but remove processed data files
        pdata1_dir = joinpath(outputexpt_str, "pdata", "1")
        if isdir(pdata1_dir)
            println("Cleaning processed data files in pdata/1...")
            processed_data_files = ["1r", "1i", "2rr", "2ri", "2ir", "2ii", "3rrr", "3iii"]
            for file in readdir(pdata1_dir)
                if file in processed_data_files
                    rm(joinpath(pdata1_dir, file), force=true)
                    println("  Removed: $file")
                end
            end
        else
            mkpath(pdata1_dir)
            println("Created pdata/1 directory")
        end
        
        # Remove other pdata subdirectories
        for dir in readdir(joinpath(outputexpt_str, "pdata"))
            if dir != "1" && isdir(joinpath(outputexpt_str, "pdata", dir))
                rm(joinpath(outputexpt_str, "pdata", dir), recursive=true)
                println("Removed pdata/$dir directory")
            end
        end
    else
        mkpath(joinpath(outputexpt_str, "pdata", "1"))
        println("Created pdata/1 directory")
    end
    
    # Function to read binary data with correct endianness and type
    function read_bruker_binary(filepath, datatype, endianness)
        open(filepath, "r") do f
            buffer = read(f)
            if endianness == "big"
                data = Array{datatype}(reinterpret(datatype, buffer))
                # Convert from big-endian to native if needed
                if ENDIAN_BOM == 0x04030201  # If system is little-endian
                    for i in 1:length(data)
                        data[i] = bswap(data[i])
                    end
                end
            else  # little-endian
                data = Array{datatype}(reinterpret(datatype, buffer))
                # Convert from little-endian to native if needed
                if ENDIAN_BOM == 0x01020304  # If system is big-endian
                    for i in 1:length(data)
                        data[i] = bswap(data[i])
                    end
                end
            end
            return data
        end
    end
    
    # Initialize arrays to store all input data and determine the minimum length
    all_data = []
    min_length = typemax(Int)
    
    # Load all input data files
    for (i, expt) in enumerate(inputexpts_str)
        datapath = joinpath(expt, data_file)
        println("Loading $(expt)/$data_file")
        
        if !isfile(datapath)
            throw(ArgumentError("Data file not found: $datapath"))
        end
        
        # Read data
        current_data = read_bruker_binary(datapath, datatype, endianness)
        push!(all_data, current_data)
        
        # Update minimum length
        min_length = min(min_length, length(current_data))
    end
    
    # Warn if truncation is needed
    if any(length(d) != min_length for d in all_data)
        println("*** WARNING! Data files have different sizes and will be truncated! ***")
        for (i, d) in enumerate(all_data)
            println("  $(inputexpts_str[i]) size: $(length(d))")
        end
        println("  Using size: $min_length")
    end
    
    # Sum in floating point, so that non-integer weights work for integer data
    summed = zeros(Float64, min_length)
    for (i, (current_data, weight)) in enumerate(zip(all_data, weights))
        summed .+= current_data[1:min_length] .* weight
        println("Applied weight $(weight) to experiment $(inputexpts_str[i])")
    end
    result_data = datatype <: Integer ? round.(datatype, summed) : summed
    
    # Write output in the format of the first experiment
    outdatapath = joinpath(outputexpt_str, data_file)
    println("Writing output to $outdatapath")
    
    open(outdatapath, "w") do f
        # Convert back to big-endian if needed
        output_data = copy(result_data)
        if endianness == "big" && ENDIAN_BOM == 0x04030201  # If system is little-endian
            for i in 1:length(output_data)
                output_data[i] = bswap(output_data[i])
            end
        elseif endianness == "little" && ENDIAN_BOM == 0x01020304  # If system is big-endian
            for i in 1:length(output_data)
                output_data[i] = bswap(output_data[i])
            end
        end
        
        write(f, reinterpret(UInt8, output_data))
    end
    
    # Calculate total scans
    total_scans = 0
    for (i, expt) in enumerate(inputexpts_str)
        acqus_path = joinpath(expt, "acqus")
        if isfile(acqus_path)
            ns_value = 0
            open(acqus_path, "r") do f
                for line in eachline(f)
                    m = match(r"##\$NS=\s*(\d+)", line)
                    if m !== nothing
                        ns_value = parse(Int, m.captures[1])
                        break
                    end
                end
            end
            weighted_scans = ns_value * abs(weights[i])
            total_scans += weighted_scans
            println("Scans from $(expt): $ns_value × $(abs(weights[i])) = $weighted_scans")
        end
    end
    
    # Update the NS value in the output acqus file
    acqufiles = ["acqus", "acqu"]
    for af in acqufiles
        output_acqus_path = joinpath(outputexpt_str, af)
        if isfile(output_acqus_path)
            acqus_content = read(output_acqus_path, String)
            updated_acqus = replace(acqus_content, r"##\$NS=\s*\d+" => "##\$NS= $(round(Int, total_scans))")
            write(output_acqus_path, updated_acqus)
        end
    end
    println("Updated number of scans to $(round(Int, total_scans))")
    
    # Create or update title file
    title_path = joinpath(outputexpt_str, "pdata", "1", "title")
    
    # Create the title content
    sum_info = "\nSum of experiments: "
    for (i, (expt, w)) in enumerate(zip(inputexpts_str, weights))
        sum_info *= "$expt (weight: $w)"
        if i < length(inputexpts_str)
            sum_info *= " + "
        end
    end
    sum_info *= "\n"
    
    # Check if title file exists, and create or append accordingly
    if isfile(title_path)
        open(title_path, "a") do f
            write(f, sum_info)
        end
    else
        mkpath(dirname(title_path))
        open(title_path, "w") do f
            write(f, sum_info)
        end
    end
    println("Updated title file with experiment summary")
    
    println("=== Operation completed successfully ===")
    return outputexpt_str
end