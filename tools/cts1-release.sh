#!/bin/bash -l

# export dry_run=1

[[ -z $package ]] && package="unknown"
[[ -z $buildflavor ]] && buildflavor="unknown"
[[ -z $rttversion ]] && rttversion="unknown"
[[ -z $access_queue ]] && access_queue="unknown"

export steps="config build test"
cmd="sbatch -J r${package:0:3}-${buildflavor: -12}-$rttversion $access_queue \
-t 4:00:00 -N 1 \
-o ${source_prefix:-notset}/logs/release-$buildflavor-$rttversion.log \
${draco_script_dir:-notset}/release.msub"
echo -e "\nConfigure, Build and Test $buildflavor-$rttversion version of $package."
echo "$cmd"

jobid=$(eval "${cmd}" &)
sleep 1m

# trim extra whitespace from number
jobid="${jobid//[^0-9]/}"
export jobids="$jobid $jobids"
