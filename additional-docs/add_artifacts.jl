# don't actually run this - it's a collection of the commands used to upload example data to github gists
# and update the Artifacts.toml file. ArtifactUtils is not a dependency of NMRTools, so
# install it into a separate environment first: `] add ArtifactUtils`, then `using ArtifactUtils`

artifact_id = artifact_from_directory("/Users/chris/Documents/Work/NMRTools/additional-data/1D_1H")
gist = upload_to_gist(artifact_id)
add_artifact!("Artifacts.toml", "1D_1H", gist, lazy=true)

artifact_id = artifact_from_directory("/Users/chris/Documents/Work/NMRTools/additional-data/1D_19F")
gist = upload_to_gist(artifact_id)
add_artifact!("Artifacts.toml", "1D_19F", gist, lazy=true)

artifact_id = artifact_from_directory("/Users/chris/Documents/Work/NMRTools/additional-data/1D_19F_titration")
gist = upload_to_gist(artifact_id)
add_artifact!("Artifacts.toml", "1D_19F_titration", gist, lazy=true)

artifact_id = artifact_from_directory("/Users/chris/Documents/Work/NMRTools/additional-data/2D_HN")
gist = upload_to_gist(artifact_id)
add_artifact!("Artifacts.toml", "2D_HN", gist, lazy=true)

artifact_id = artifact_from_directory("/Users/chris/Documents/Work/NMRTools/additional-data/2D_HN_titration")
gist = upload_to_gist(artifact_id)
add_artifact!("Artifacts.toml", "2D_HN_titration", gist, lazy=true)

artifact_id = artifact_from_directory("/Users/chris/Documents/Work/NMRTools/additional-data/pseudo2D_T2")
gist = upload_to_gist(artifact_id)
add_artifact!("Artifacts.toml", "pseudo2D_T2", gist, lazy=true)

artifact_id = artifact_from_directory("/Users/chris/Documents/Work/NMRTools/additional-data/pseudo2D_XSTE")
gist = upload_to_gist(artifact_id)
add_artifact!("Artifacts.toml", "pseudo2D_XSTE", gist, lazy=true)

artifact_id = artifact_from_directory("/Users/chris/Documents/Work/NMRTools/additional-data/pseudo3D_HN_R2")
gist = upload_to_gist(artifact_id)
add_artifact!("Artifacts.toml", "pseudo3D_HN_R2", gist, lazy=true)

artifact_id = artifact_from_directory("/Users/chris/Documents/Work/NMRTools/additional-data/pseudo3D_kinetics")
gist = upload_to_gist(artifact_id)
add_artifact!("Artifacts.toml", "pseudo3D_kinetics", gist, lazy=true)

artifact_id = artifact_from_directory("/Users/chris/Documents/Work/NMRTools/additional-data/3D_HNCA")
gist = upload_to_gist(artifact_id)
add_artifact!("Artifacts.toml", "3D_HNCA", gist, lazy=true)

artifact_id = artifact_from_directory("/Users/chris/Documents/Work/NMRTools/additional-data/3D_HNCA_pipe")
gist = upload_to_gist(artifact_id)
add_artifact!("Artifacts.toml", "3D_HNCA_pipe", gist, lazy=true)

artifact_id = artifact_from_directory("/Users/chris/Documents/Work/NMRTools/additional-data/3D_HNCA_ucsf")
gist = upload_to_gist(artifact_id)
add_artifact!("Artifacts.toml", "3D_HNCA_ucsf", gist, lazy=true)

