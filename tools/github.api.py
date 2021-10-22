#!/usr/bin/env python
# ------------------------------------------------------------------------------------------------ #
# File: github.api.py
# Note: Copyright 2020 Triad National Security, LLC., All rights reserved.
#
# Sample use:
# ---------------------------------------
# python ./github.api.py --user <github username> \
#        --token <api token key> <date of last release, YYYYMMDD>
# python ./github.api.py --user KineticTheory \
#        --token abcdef1234567890abcdef1234567890abcdef12 20190713
# ------------------------------------------------------------------------------------------------ #

import requests
# import json
import re
# import time
from datetime import date
import argparse

# Defaults
start_date = 20000101  # YYYYMMDD

# First (only) arg is the cut-off date.  Only more recently merged PRs will be printed.
parser = argparse.ArgumentParser()
parser.add_argument("date_of_last_release", help="Date of last release (YYYYMMDD)")
parser.add_argument("--user", help="github user name")
parser.add_argument("--token",
                    help="github user token (from github user profile developer settings)")
args = parser.parse_args()
start_date = int(args.date_of_last_release)
username = args.user
git_token = args.token

if username is None:
    print("\nERROR: You must provide values for --user, --token, and the date of the last release.")
    print("         No value for --user was provided.\n")
    parser.print_help()
    quit()

if git_token is None:
    print("\nERROR: You must provide values for --user, --token, and the date of the last release.")
    print("         No value for --token was provided.\n")
    parser.print_help()
    quit()

# This API url gets all of the closed pull resquests in draco, this generally looks like
# https://api.github.com/DATA/ORG/REPO_NAME/DESCRIPTOR, see the github API for more information
api_url = 'https://api.github.com/repos/lanl/Draco/pulls?state=closed'
result = requests.get(api_url, auth=(username, git_token))
result_json = result.json()

# input the start date as YYYYMMDD for ease in comparing
lryear = int(start_date / 10000)
lrmonth = int((start_date - lryear * 10000) / 100)
lrday = int((start_date - lryear * 10000 - lrmonth * 100))
lastreleasedate = date(lryear, lrmonth, lrday)
start_date_str = lastreleasedate.strftime("%B %d %Y")

re_date_str = re.compile("([0-9]{4})-([0-9]{2})-([0-9]{2})")
all_done = False

# ------------------------------------------------------------------------------------------------ #
# Pull Requests
# ------------------------------------------------------------------------------------------------ #

# parse the first json result page, onlt print PR data for PRs that were merged
print("\nShowing PRs merged to develop dated after {0}:\n".format(start_date_str))

for entry in result_json:
    if entry == "message":
        if result_json["message"] == "Bad credentials":
            raise Exception('Bad credentials. Check GitHub moniker and token string.')
    if entry["merged_at"]:
        merge_date = entry["merged_at"]
        base = entry["base"]
        # restrict report to PRs that are:
        # - merged (entry 'merge_date' exists).
        # - merged to develop. (not 100% about the validity of this check)
        if (merge_date and base and base["label"] == "lanl:develop"):
            str_result = re_date_str.findall(merge_date)[0]
            number_date = 10000 * int(str_result[0]) + 100 * int(str_result[1]) + int(str_result[2])
            if (number_date > start_date):
                print("* [PR #{0} {1}]({2})".format(entry["number"], entry["title"],
                                                    entry["html_url"]))
            else:
                all_done = True
                break

# if there are more results than will fit in the first page of results, then continue
while 'next' in result.links.keys() and not all_done:
    result = requests.get(result.links['next']['url'], auth=(username, git_token))
    result_json = result.json()
    for entry in result_json:
        merge_date = entry["merged_at"]
        base = entry["base"]
        # restrict report to PRs that are:
        # - merged (entry 'merge_date' exists).
        # - merged to develop. (not 100% about the validity of this check)
        if (merge_date and base and base["label"] == "lanl:develop"):
            str_result = re_date_str.findall(merge_date)[0]
            number_date = 10000 * int(str_result[0]) + 100 * int(str_result[1]) + int(str_result[2])
            if (number_date > start_date):
                print("* [PR #{0} {1}]({2})".format(entry["number"], entry["title"],
                                                    entry["html_url"]))
            else:
                all_done = True
                break

# ------------------------------------------------------------------------------------------------ #
# Issues - Fixed
# ------------------------------------------------------------------------------------------------ #

# this API url gets all of the closed pull resquests in draco, this generally looks like
# https://api.github.com/DATA/ORG/REPO_NAME/DESCRIPTOR, see the github API for more information

api_url = 'https://api.github.com/repos/lanl/Draco/issues?state=closed'
result = requests.get(api_url, auth=(username, git_token))
result_json = result.json()
all_done = False

# parse the first json result page, onlt print issue data for issues that were closed
print("\nShowing issues closed after {0}:\n".format(start_date_str))
for entry in result_json:
    if "pull_request" not in entry:  # this entry only exists for PRs not issues
        closed_date = entry["closed_at"]
        if (closed_date):
            str_result = re_date_str.findall(closed_date)[0]
            number_date = 10000 * int(str_result[0]) + 100 * int(str_result[1]) + int(str_result[2])
            if (number_date > start_date):
                labels = entry["labels"]
                if(labels):
                    print("* [Github {0} #{1} {2}]({3})".format(labels[0]["name"], entry["number"],
                                                                entry["title"], entry["html_url"]))
                else:
                    print("* [Github issue #{0} {1}]({2})".format(entry["number"], entry["title"],
                                                                  entry["html_url"]))
            else:
                all_done = True
                break

# The rest of this is broken.  I need to spend more time with it but not right now.

# if(not all_done):
#     # if there are more results than will fit in the first page of results, then
#     # continue
#     while 'next' in result.links.keys():
#         result = requests.get(result.links['next']['url'], auth=(username, git_token))
#         result_json = result.json()
#         for entry in result_json:
#             print(entry)
#             if not "pull_request" in entry: # this entry only exists for PRs not issues
#                 closed_date = entry["cloased_at"]
#                 if (closed_date):
#                     str_result = re_date_str.findall(closed_date)[0]
#                     number_date = 10000 * int(str_result[0]) + 100 * int(str_result[1]) \
#                                   + int(str_result[2])
#                     if (number_date > start_date):
#                         print("\"Github issue #{0} {1}\":{2}".format(entry["iid"],
#                                                                      entry["title"],
#                                                                      entry["web_url"]))
#                     else:
#                         break

# ------------------------------------------------------------------------------------------------ #
# Issues - Still Open (label:bug)
# ------------------------------------------------------------------------------------------------ #

# this API url gets all of the closed pull resquests in draco, this generally looks like
# https://api.github.com/DATA/ORG/REPO_NAME/DESCRIPTOR, see the github API for more information

api_url = 'https://api.github.com/repos/lanl/Draco/issues?state=open'
result = requests.get(api_url, auth=(username, git_token))
result_json = result.json()
all_done = False

# parse the first json result page, onlt print issue data for issues that were closed
print("\nShowing open issues (most recently reported):\n")
for entry in result_json:
    if "pull_request" not in entry:  # this entry only exists for PRs not issues
        # only print issues labeled as "bug"
        if entry["labels"]:
            for label in entry["labels"]:
                if label["name"] == "bug":
                    print("* [Github bug #{0} {1}]({2})".format(entry["number"], entry["title"],
                                                                entry["html_url"]))

# ------------------------------------------------------------------------------------------------ #
# End github.api.py
# ------------------------------------------------------------------------------------------------ #
