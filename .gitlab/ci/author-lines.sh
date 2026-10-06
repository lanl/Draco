#!/usr/bin/env bash
# Count surviving selected-file lines by author at HEAD (including comments/blanks).
# Usage: author-lines.sh -s <repository-root>
# Output: author-lines.tsv in the directory where you invoke this script.
set -euo pipefail

# Global inclusion globs. Quote entries to prevent shell expansion here.
FILE_PATTERNS=(
    '*.c' '*.cc' '*.cpp' '*.cxx' '*.h' '*.hh' '*.hpp' '*.hxx'
    '*.in' '*.json' '*.md' '*.txt'
)

# Bash extended regular expressions matched against repository-relative paths.
# A match to ANY entry excludes the file. Empty means no exclusions.
# Uncomment or add entries as needed (regex syntax, not shell glob syntax).
EXCLUDE_REGEXES=(
    '(^|/)[^/]*terminal[^/]*/'
    # '(^|/)(third_party|vendor|generated)/'
    # '^docs/'
    # '(^|/)CMakeLists\.txt$'
)

# One entry per person:
# 'Canonical name|canonical@email|alias@email|another-alias@email'
# The canonical email is also registered, so different names using it merge.
# Use bare emails (no angle brackets). Matching is exact and case-sensitive.
# Entries cannot contain literal pipes, tabs, or newlines in names/emails.
# Git applies any repository mailmap first; register the emails blame emits.
AUTHOR_REGISTRY=(
  'Alex Long|along@lanl.gov|arlong.ne@gmail.com|along@darwin-fe1.lanl.gov|along@9c8f39b7-3d63-494d-b0b1-10263151bafb'
  'Alexis Maldonado|amaldonado@lanl.gov'
  'Allan Wollaber|wollaber@2e4270bf-7bdd-41b3-a9ef-8cbd4d332402|wollaber@9c8f39b7-3d63-494d-b0b1-10263151bafb'
  'Andrew Thomas Till|till@lanl.gov|attom@users.noreply.github.com'
  'Ben Ryan|brryan@lanl.gov'
  'Daniel Holladay|dholladay00@lanl.gov'
  'David Dixon|david8dixon@gmail.com'
  'Gabe Rockefeller|gaber@lanl.gov|gaber@2e4270bf-7bdd-41b3-a9ef-8cbd4d332402|gaber@9c8f39b7-3d63-494d-b0b1-10263151bafb'
  'HyeongKae (Ryosuke) Park|hkpark@lanl.gov|hkpark@9c8f39b7-3d63-494d-b0b1-10263151bafb|hkpark@ml-fey1.lanl.gov|hkpark@pi-fey.lanl.gov'
  'Jae Chang|jhchang@lanl.gov|jhchang@2e4270bf-7bdd-41b3-a9ef-8cbd4d332402'
  'James S. Warsa|warsa@lanl.gov|warsa@2e4270bf-7bdd-41b3-a9ef-8cbd4d332402'
  'Jeff Densmore|jdd@9c8f39b7-3d63-494d-b0b1-10263151bafb'
  'Jeffrey Robert Haack|haack@lanl.gov'
  'Joseph Michael Coale|jmcoale@lanl.gov'
  'Hans Ruediger Hammer|hrhammer@lanl.gov'
  'Irina Sagert|sagert@lanl.gov|sagert@darwin-fe1.lanl.gov'
  'Kelly (KT) Thompson|kgt@lanl.gov|kellyt@2e4270bf-7bdd-41b3-a9ef-8cbd4d332402|107638@2e4270bf-7bdd-41b3-a9ef-8cbd4d332402|KineticTheory@users.noreply.github.com|kellyt@9c8f39b7-3d63-494d-b0b1-10263151bafb'
  'Kendra Phyllis Long|keadyk@lanl.gov|keadyk@2e4270bf-7bdd-41b3-a9ef-8cbd4d332402|long60@rzadams1017.llnl.gov|long60@rzansel61.coral.llnl.gov|keadyk@darwin-fe1.lanl.gov|keadyk@9c8f39b7-3d63-494d-b0b1-10263151bafb'
  'Kent G. Budge|kgbudge@lanl.gov|kgbudge@2e4270bf-7bdd-41b3-a9ef-8cbd4d332402'
  'Katherine J. Wang|kwang@2e4270bf-7bdd-41b3-a9ef-8cbd4d332402'
  'Kyle Stuart Beling|kbeling@lanl.gov'
  'Mathew Cleveland|clevelam@lanl.gov|cleveland@lanl.gov|clevelam@2e4270bf-7bdd-41b3-a9ef-8cbd4d332402|clevelam@9c8f39b7-3d63-494d-b0b1-10263151bafb'
  'Peter Ahrens|pahrens@2e4270bf-7bdd-41b3-a9ef-8cbd4d332402'
  'Paul Talbot|talbotp@2e4270bf-7bdd-41b3-a9ef-8cbd4d332402'
  'Rob Lowrie|lowrie@lanl.gov|lowrie@pn1244044.lanl.gov|lowrie@9c8f39b7-3d63-494d-b0b1-10263151bafb'
  'Ryan Thomas Wollaeger|wollaeger@lanl.gov|wollaege@9c8f39b7-3d63-494d-b0b1-10263151bafb|wollaege@darwin-fe2.lanl.gov'
  'Timothy Kelley|tkelley@lanl.gov|timmah@users.noreply.github.com'
  'Todd Urbatsch|tmonster@lanl.gov|tmonster@9c8f39b7-3d63-494d-b0b1-10263151bafb'
)

usage() {
    printf 'Usage: %s -s <repository-root>\n' "${0##*/}"
}

repo=
while getopts ':s:h' opt; do
    case "$opt" in
        s) repo=$OPTARG ;;
        h) usage; exit 0 ;;
        :) printf 'Option -%s requires an argument.\n' "$OPTARG" >&2
           usage >&2; exit 2 ;;
        \?) printf 'Unknown option: -%s\n' "$OPTARG" >&2
            usage >&2; exit 2 ;;
    esac
done
shift "$((OPTIND - 1))"
if [[ -z $repo ]] || (( $# != 0 )); then
    usage >&2
    exit 2
fi

# Reject invalid exclusion regexes instead of silently ignoring them.
for regex in "${EXCLUDE_REGEXES[@]}"; do
    if [[ '' =~ $regex ]]; then
        :
    else
        status=$?
        if (( status == 2 )); then
            printf 'Invalid exclusion regex: %s\n' "$regex" >&2
            exit 2
        fi
    fi
done

# Capture the caller's working directory before changing to the repository.
output="$PWD/author-lines.tsv"
root=$(git -C "$repo" rev-parse --show-toplevel)
cd "$root"
revision=$(git rev-parse --verify 'HEAD^{commit}')

# Publish only after the entire pipeline succeeds.
tmp=$(mktemp -- "${output}.tmp.XXXXXX")
trap 'rm -f -- "$tmp"' EXIT

git ls-tree -r -z --name-only "$revision" |
while IFS= read -r -d '' file; do
    included=false
    for pattern in "${FILE_PATTERNS[@]}"; do
        # Unquoted RHS intentionally enables glob matching inside [[ ... ]].
        if [[ $file == "$pattern" ]]; then
            included=true
            break
        fi
    done
    [[ $included == true ]] || continue

    for regex in "${EXCLUDE_REGEXES[@]}"; do
        # Unquoted RHS intentionally enables regex matching.
        if [[ $file =~ $regex ]]; then
            continue 2
        fi
    done

    # -w ignores whitespace changes; -M detects moves within a file.
    git blame --line-porcelain -w -M "$revision" -- "$file" || exit 1
done |
AUTHOR_REGISTRY_DATA="$(printf '%s\n' "${AUTHOR_REGISTRY[@]}")" awk '
    BEGIN {
        n = split(ENVIRON["AUTHOR_REGISTRY_DATA"], entries, "\n")
        for (i = 1; i <= n; i++) {
            if (entries[i] == "") continue
            fields = split(entries[i], parts, "[|]")
            if (fields < 2 || parts[1] == "" || parts[2] == "" ||
                entries[i] ~ /\t/) {
                print "Invalid AUTHOR_REGISTRY entry: " entries[i] > "/dev/stderr"
                exit 2
            }
            canonical = parts[1] "\t<" parts[2] ">"
            for (j = 2; j <= fields; j++) {
                email = parts[j]
                if (email == "" || email ~ /[<>]/ ||
                    (email in identities && identities[email] != canonical)) {
                    print "Invalid or conflicting registry email: " email > "/dev/stderr"
                    exit 2
                }
                identities[email] = canonical
            }
        }
    }
    /^author / { name = substr($0, 8) }
    /^author-mail / {
        original_email = substr($0, 13)
        email = original_email
        sub(/^</, "", email)
        sub(/>$/, "", email)
        identity = (email in identities) ? identities[email] : name "\t" original_email
        counts[identity]++
    }
    END {
        for (author in counts)
            printf "%d\t%s\n", counts[author], author
    }
' |
LC_ALL=C sort -t $'\t' -k1,1nr -k2,2 -k3,3 |
{
    printf 'Lines\tAuthor\temail\n'
    cat
} > "$tmp"

mv -f -- "$tmp" "$output"
printf 'Wrote %s (revision %s)\n' "$output" "$revision" >&2
