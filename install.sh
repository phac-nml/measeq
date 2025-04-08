#!/bin/bash
set -o pipefail

### Set Variables ###
#####################
SCRIPT_DIR=$(dirname "$(realpath "$0")")
prefix=""
name=""
CONDA_PATH=""
CONDA_NAME=""
prefix_chosen=false
tool="conda"
HELP="
USAGE:\n
./install.sh \n\n

Optional:\n
    -p | --prefix <DIR>         : Prefix to where the conda environment will be created.\n
    -n | --name <NAME>          : Name of the conda environment to be created.\n
    --use-mamba                 : Use mamba instead of conda for the installation.\n\n
"

# Check if Conda exists
if ! command -v conda &> /dev/null; then
    echo "ERROR: Conda is not installed. Exiting."
    exit 1
fi

# Initialize conda
eval "$(conda shell.bash hook)"

### Arg Parsing and Validation ###
##################################
while [[ $# -gt 0 ]]; do
    case $1 in
        -p|--prefix)
            prefix=$(realpath "$2")
            shift 2
            ;;
        -n|--name)
            name="$2"
            shift 2
            ;;
        --use-mamba)
            tool="mamba"
            shift 1
            ;;
        -h|--help)
            echo -e $HELP
            exit 0
            ;;
        *)
            echo "ERROR: Unknown option detected ($1)"
            exit 1
            ;;
    esac
done

### Install from environment file ###
#####################################
cd $SCRIPT_DIR

# Check what arguments were provided
if [[ -n "$prefix" ]] && [[ -n "$name" ]]; then
    echo "ERROR: Please provide only a prefix or a name for the conda environment, not both."
    exit 1
elif [[ -n "$prefix" ]] && [[ -z "$name" ]]; then
    CONDA_PATH="$prefix"
elif [[ -z "$prefix" ]] && [[ -n "$name" ]]; then
    CONDA_NAME="$name"
else
    CONDA_NAME="measeq"
fi

# Run installation
if [[ ! -z "$CONDA_PATH" ]]; then
    $tool env create -f "$SCRIPT_DIR/measeq/environment.yml" -p "$CONDA_PATH"
    prefix_chosen=true
else 
    $tool env create -f "$SCRIPT_DIR/measeq/environment.yml" -n "$CONDA_NAME"
    CONDA_PATH=$(conda env list | grep "^$CONDA_NAME" | awk '{print $2}')
fi

### Post installation file organization ###
###########################################
# Move reference data and scripts
cp -r $SCRIPT_DIR/measeq/src/measeq/ "$CONDA_PATH/lib"

# Move scripts
chmod +x $SCRIPT_DIR/measeq/measeq
cp $SCRIPT_DIR/measeq/measeq "$CONDA_PATH/bin"

### Test Installation ###
#########################
conda activate "$CONDA_PATH"
measeq -h > /dev/null 2>&1
if [[ $? -eq 0 ]]; then
    echo -e "\n\nINFO: Conda environment created and is available using the following command:\n\n
    \e[32mconda activate $( [ "$prefix_chosen" = true ] && echo "$CONDA_PATH" || echo "$CONDA_NAME" )\e[0m\n\n"
else
    echo "ERROR: There was an error with the installation."
    exit 1
fi
conda deactivate
