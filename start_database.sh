#!/usr/bin/env bash

# Define the path to your Docker Compose YAML file
COMPOSE_FILE="mongodb/docker-compose.yml"

# Check if Docker Compose is installed
if ! command -v docker-compose &>/dev/null; then
    echo "Docker Compose is not installed. Please install it first."
    exit 1
fi

# Navigate to the directory containing the Docker Compose file
cd "$(dirname "$COMPOSE_FILE")" || exit 1

# Start the Docker containers using Docker Compose
docker-compose up -d

# Check if the containers were started successfully
if [ $? -ne 0 ]; then
    echo "Failed to start Docker containers."
    exit 1
fi

# Exit the script
exit 0