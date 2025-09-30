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



