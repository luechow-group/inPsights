
function(add_to_source_group name_prefix folder_location variable_name)
    foreach(child ${${variable_name}})
        get_filename_component(directory_name ${child} DIRECTORY)
        string(REPLACE "${folder_location}" "${name_prefix}" prefixed_groupname ${directory_name})
        string(REPLACE "." "\\" groupname ${prefixed_groupname})
        source_group("${groupname}" FILES ${child})
    endforeach()
endfunction()

