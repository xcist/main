# Copyright 2024, GE Precision HealthCare. All rights reserved. See https://github.com/xcist/main/tree/master/license

"""
Aim
    Checks whether file 1 is newer than file 2

Inputs
    FirstFile   A string containing the first filename
    SecondFile  A string containing the second filename

Outputs
    ReturnValue 0 if FirstFile does not exist
                1 if FirstFile exists and SecondFile does not exist
                  or if FirstFile is newer than SecondFile
"""
import os

def first_file_newer(first_file, second_file):
    # If first_file does not exist, return False. If first_file exists and second_file does not,
    # return True. Otherwise, return True if first_file is newer than second_file.

    if not os.path.exists(first_file):
        return False

    if os.path.exists(first_file) and not os.path.exists(second_file):
        return True

    res1 = os.path.getmtime(first_file)
    res2 = os.path.getmtime(second_file)

    return first_date_is_more_recent(res1, res2)

def first_date_is_more_recent(date1, date2):
    t1 = monotonic_with_time(date1)
    t2 = monotonic_with_time(date2)

    return t1 > t2

def monotonic_with_time(date):
    date = date.replace(':', ' ').replace('-', ' ')
    split = date.split()

    if len(split) != 6:
        raise ValueError(f'Incompatible date format: {date}')

    out = int(split[2]) * 12 * 31 * 24 * 3600  # years --> seconds
    out += (numeric_month(split[1]) - 1) * 31 * 24 * 3600  # months --> seconds
    out += (int(split[0]) - 1) * 24 * 3600  # days --> seconds
    out += int(split[3]) * 3600  # hours --> seconds
    out += int(split[4]) * 60  # minutes --> seconds
    out += int(split[5])  # seconds

    return out

def numeric_month(month):
    months = {
        'Jan': 1, 'Feb': 2, 'Mar': 3, 'Apr': 4, 'May': 5, 'Jun': 6,
        'Jul': 7, 'Aug': 8, 'Sep': 9, 'Oct': 10, 'Nov': 11, 'Dec': 12
    }
    return months[month]