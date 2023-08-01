#!/bin/bash

# Check if command line argument is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 DIRECTORY_TO_MONITOR"
    exit 1
fi

# The directory to monitor
DIRECTORY_TO_MONITOR=$1

echo "Monitoring $DIRECTORY_TO_MONITOR"

# Inotifywait options: monitor for create, modify events
INOTIFY_EVENTS="create,modify"

# Run indefinitely
while true
do
    # Use inotifywait to wait for changes in DIRECTORY_TO_MONITOR and its subdirectories, outputting only the event and filename
    change=$(inotifywait -rq -e "$INOTIFY_EVENTS" --format '%w%f %e' "$DIRECTORY_TO_MONITOR")

    # Parse the output into the file path and event
    filepath=$(echo $change | awk '{print $1}')
    filepath="${filepath%.*}"
    event=$(echo $change | awk '{print $2}')

    echo "candidate $filepath"

    # Check if the changed file is a .ui file
    if [[ $filepath == *.ui ]]
    then
        # If a .ui file was changed, run pyuic5 on it
        echo "File $filepath has changed, running pyuic5..."

        # Get the directory and base filename without the extension
        dir=$(dirname "$filepath")
        base=$(basename "$filepath" .ui)

        # Run pyuic5 and output to a .py file with the "Ui_" prefix
        pyuic5 "$filepath" -o "$dir/Ui_${base}.py"
    fi
done
