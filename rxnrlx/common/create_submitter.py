

def create_submit_script(account_info:dict, command:str):
    """
    Create submit script from user-specified inputs
    """


    submit_script = ["#!/bin/bash -l", ""]

    # iterate through the specified parameters and build the script
    for key, val in account_info:
        submit_script.append(f"#SBATCH --{key}={val}")
    
    # append final command
    submit_script.append("")
    submit_script.append(command)

    with open("submit.script", "w") as f: 
        f.write("\n".join(submit_script))
