#!/usr/bin/env bash

venv_path="./venv"

if [ -d "$venv_path" ]; then
    source "$venv_path/bin/activate"
    echo "Virtual environment activated"
    python3.9 pathogen_retrieve/pathogen_data_mining.py
    python3.9 visualize_data/plot_distribution.py
    deactivate
else
    echo "No venv directory in: $venv_path"
    exit 1
fi
