#!/usr/bin/env python
# ------------------------------------------------------------------------------------------------ #
# File: gitlab.api.py
# Note: Copyright (C) 2022 Triad National Security, LLC., All rights reserved.
#
# Sample use:
# ---------------------------------------
# python ./gitlab.api.py --user <gitlab username> \
#        --token <api token key> <date of last release, YYYYMMDD>
# python ./gitlab.api.py --user kellyt \
#        --token abcdef1234567890abcdef1234567890abcdef12 20190713
# ------------------------------------------------------------------------------------------------ #

import requests
import re
import os
# import time
from datetime import date
import argparse

# Defaults

start_date = 20000101  # YYYYMMDD

# First (only) arg is the cut-off date.  Only more recently merged PRs will be printed.
parser = argparse.ArgumentParser()
parser.add_argument("date_of_last_release", help="Date of last release (YYYYMMDD)")
parser.add_argument("--user", help="gitlab user name")
parser.add_argument("--token",
                    help="gitlab user token (from gitlab user profile developer settings)")
args = parser.parse_args()
start_date = int(args.date_of_last_release)
username = args.user
git_token = args.token

if username is None:
    print("\nERROR: You must provide values for --user, --token, and the date of the last release.")
    print("         No value for --user was provided.\n")
    parser.print_help()
    quit()

if (username is None or git_token is None):
    print("\nERROR: You must provide values for --user, --token, and the date of the last release.")
    print("         No value for --token was provided.\n")
    parser.print_help()
    quit()

os.environ['no_proxy'] = '*'

# For the gitlab API you have to use the project ID instead of the project name I couldn't find this
# anywhere on our gitlab page but it shows up if you list all the projects through the API
# MAC: I found this by looking at the HTML source of the draco gitlab page (see
# https://stackoverflow.com/questions/39559689/where-do-i-find-the-project-id-for-the-gitlab-api)
# for instructions.
draco_id = "503"  # Draco
api_url = 'https://re-git.lanl.gov/api/v4/projects/{0}/merge_requests?'.format(draco_id)
api_url += 'state=merged&private_token={0}'.format(git_token)

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
# Merge Requests
# ------------------------------------------------------------------------------------------------ #

# parse the first json result page, onlt print PR data for PRs that were merged
print("\nShowing MRs merged to develop dated after {0}:\n".format(start_date_str))
for entry in result_json:
    merge_date = entry["merged_at"]
    # restrict report to PRs that are:
    # - merged (entry 'merge_date' exists).
    # - merged to develop. (not 100% about the validity of this check)
    if (merge_date and entry["target_branch"] == "develop"):
        str_result = re_date_str.findall(merge_date)[0]
        number_date = 10000 * int(str_result[0]) + 100 * int(str_result[1]) + int(str_result[2])
        if (number_date > start_date):
            print("* [MR !{0} {1}]({2})".format(entry["iid"], entry["title"],
                                                entry["web_url"]))
        else:
            all_done = True
            break

# if there are more results than will fit in the first page of results, then continue
while 'next' in result.links.keys() and not all_done:
    result = requests.get(result.links['next']['url'], auth=(username, git_token))
    result_json = result.json()
    for entry in result_json:
        merge_date = entry["merged_at"]
        # restrict report to PRs that are:
        # - merged (entry 'merge_date' exists).
        # - merged to develop. (not 100% about the validity of this check)
        if (merge_date and entry["target_branch"] == "develop"):
            str_result = re_date_str.findall(merge_date)[0]
            number_date = 10000 * int(str_result[0]) + 100 * int(str_result[1]) + int(str_result[2])
            if (number_date > start_date):
                print("* [MR !{0} {1}]({2})".format(entry["iid"], entry["title"],
                                                    entry["web_url"]))
            else:
                all_done = True
                break

# ------------------------------------------------------------------------------------------------ #
# Issues - Fixed
# ------------------------------------------------------------------------------------------------ #

# this API url gets all of the closed merge resquests in Draco, this generally looks like
# https://api.github.com/DATA/ORG/REPO_NAME/DESCRIPTOR, see the github API for more information

api_url = 'https://re-git.lanl.gov/api/v4/projects/{0}/issues?'.format(draco_id)
api_url += 'state=closed&private_token={0}'.format(git_token)
result = requests.get(api_url, auth=(username, git_token))
result_json = result.json()
all_done = False

# parse the first json result page, onlt print issue data for issues that were closed
print("\nShowing issues closed after {0}:\n".format(start_date_str))
for entry in result_json:
    # if "pull_request" not in entry:  # this entry only exists for PRs not issues
    closed_date = entry["closed_at"]
    if (closed_date):
        str_result = re_date_str.findall(closed_date)[0]
        number_date = 10000 * int(str_result[0]) + 100 * int(str_result[1]) + int(str_result[2])
        if (number_date > start_date):
            labels = entry["labels"]
            issuetype = 'issue'
            if 'bug' in labels:
                issuetype = 'bug'
            print("* [Gitlab {0} #{1} {2}]({3})".format(issuetype, entry["iid"], entry["title"],
                                                        entry["web_url"]))
        else:
            all_done = True
            break

# ------------------------------------------------------------------------------------------------ #
# Issues - Still Open (label:bug)
# ------------------------------------------------------------------------------------------------ #

# this API url gets all of the closed pull resquests in draco, this generally looks like
# https://api.github.com/DATA/ORG/REPO_NAME/DESCRIPTOR, see the github API for more information

api_url = 'https://re-git.lanl.gov/api/v4/projects/{0}/issues?'.format(draco_id)
api_url += 'state=opened&private_token={0}'.format(git_token)
result = requests.get(api_url, auth=(username, git_token))
result_json = result.json()
all_done = False

# parse the first json result page, onlt print issue data for issues that were closed
print("\nShowing open issues (most recently reported):\n")
for entry in result_json:
    # only print issues labeled as "bug"
    labels = entry["labels"]
    if 'bug' in labels:
        print("* [Gitlab bug #{0} {1}]({2})".format(entry["iid"], entry["title"],
                                                    entry["web_url"]))

# ------------------------------------------------------------------------------------------------ #
# End
# ------------------------------------------------------------------------------------------------ #
