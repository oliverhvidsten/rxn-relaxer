""" Utility functions that are broadly applicable """

def sec_to_str(seconds):
    """Converts a numeric time value into a human readable hours:minutes:seconds format"""

    # get days and remove from total
    days = seconds // (3600 * 24)
    seconds -= (days * 3600 * 24)

    # get hours and remove from total
    hr = seconds // 3600
    seconds -= (hr * 3600)
    
    # get minutes and remove from total (only seconds should be left)
    min = seconds // 60
    seconds -= (min * 60)


    # return formatted string
    return f"{days}-{hr}:{min}:{seconds:.2f}"



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