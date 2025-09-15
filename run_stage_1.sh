#!/usr/bin/env bash  
set -euo pipefail    

GITHUB_REMOTE_NAME="${GITHUB_REMOTE_NAME:-origin}"
LOCAL_BRANCH="${LOCAL_BRANCH:-$(git rev-parse --abbrev-ref HEAD)}"
REMOTE_HOST="${REMOTE_HOST:?e.g. user@workstation}"
REMOTE_REPO_DIR="${REMOTE_REPO_DIR:?e.g. /home/user/myrepo}"
OUT_ROOT="${OUT_ROOT:-output}"        # remote output root
BUILD_CMD="${BUILD_CMD:-make -j}"     # remote build cmd #######################
RUN_CMD="${RUN_CMD:-./my_program}"    # remote run cmd

# 1) Local: Check whether everything is committed 
if [[git status --porcelain]]; then
  echo "You have uncommitted changes!" >&2
  exit 1
fi

# 2) Commit and push to GitHub
COMMIT="$(git rev-parse HEAD)" #################################################
git push "$GITHUB_REMOTE_NAME" "$LOCAL_BRANCH" --follow-tags

# 3) Verify GitHub has THIS commit
git fetch -q "$GITHUB_REMOTE_NAME"
if [[ "$(git merge-base --is-ancestor "$COMMIT" "$GITHUB_REMOTE_NAME/$LOCAL_BRANCH"; echo $?)" != "0" ]]; then
  echo "ERROR: Commit $COMMIT is not on $GITHUB_REMOTE_NAME/$LOCAL_BRANCH after push." >&2
  exit 1
fi

# 4) Trigger the workstation: build exactly this commit
ssh -o BatchMode=yes "$REMOTE_HOST" bash -s <<EOF
set -euo pipefail
cd "$REMOTE_REPO_DIR"
export GIT_COMMIT="$COMMIT"
export OUT_ROOT="$OUT_ROOT"
export BUILD_CMD='$BUILD_CMD'
export RUN_CMD='$RUN_CMD'
bash run_strict_by_commit.sh
EOF
