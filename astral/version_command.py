import subprocess
import re

def get_astral_version():
    try:
        #run the `astral --help` command
        result = subprocess.run(['astral', '--help'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        
        #save the output
        output = result.stdout + result.stderr
        
        # Define the regex pattern to match the version string from the help message
        version_pattern = re.compile(r'This is ASTRAL version (\d+\.\d+\.\d+)')
        
        #search for the version pattern in the output
        match = version_pattern.search(output)
        
        if match:
            # Extract and return the version string
            return match.group(1)
        else:
            print("Version information not found in `astral --help` output.")
            return None
    except Exception as e:
        print(f"An error occurred: {e}")
        return None

# Call the function and print the version
version = get_astral_version()
if version:
    print(version)
