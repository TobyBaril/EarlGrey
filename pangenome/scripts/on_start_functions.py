import os
import sys
import shutil
import subprocess
import urllib.request
from pathlib import Path

def running_tea(stage="Starting Earl Grey"):
    """Display ASCII art tea cup with stage name"""
    tea_art = rf"""    
          )  (
         (   ) )
         ) ( (
       _______)_
    .-'---------|  
   ( C|/\/\/\/\/|
    '-./\/\/\/\/|
      '_________'
       '-------'
    <<< {stage} >>>"""
    print(tea_art)

def convert_seconds(seconds):
    """Convert seconds to HH:MM:SS.ss format"""
    seconds = int(seconds)
    h = seconds // 3600
    m = seconds % 3600 // 60
    s = seconds % 60
    return f"{h:02d}:{m:02d}:{s:05.2f}"

def validate_parameters(config):
    """Validate and set default parameters"""
    required_params = ['genome', 'species', 'output_dir']
    
    # Check required parameters
    for param in required_params:
        if not config.get(param):
            print(f"ERROR: Required parameter '{param}' not specified")
            sys.exit(1)
    
    # Set defaults and display messages
    defaults = {
        'ProcNum': (1, "Cores Will Be Used"),
        'num': (10, "De Novo Sequences Will Be Extended Through a Maximum of {} Iterations"),
        'no_seq': (20, "{} sequences will be used in BEAT consensus generation"),
        'cluster': ('no', None),
        'softMask': ('no', None),
        'margin': ('no', None),
        'Flank': (1000, "Blast, Extend, Align, Trim Process Will Add {}bp to Each End in Each Iteration"),
        'min_seq': (3, "Blast, Extend, Align, Trim Process Will Require {} Sequences to Generate a New Consensus Sequence"),
        'heli': ('no', None)
    }
    
    for param, (default_val, message_template) in defaults.items():
        if not config.get(param):
            config[param] = default_val
            if message_template:
                print(message_template.format(default_val))
        else:
            if message_template:
                print(message_template.format(config[param]))
    
    # Handle special cases with custom messages
    if not config.get('RepSpec') and not config.get('startCust'):
        print("RepeatMasker species not specified, running Earl Grey without an initial mask with known repeats")
    else:
        print("Running Earl Grey with an initial mask with known repeats")
    
    # Clustering validation
    if config['cluster'] == 'yes':
        print("TE library consensus sequences will be clustered, be aware of the impact this can have on subfamilies and chimeric TEs!")
    else:
        print("TE library consensus sequences will not be clustered")
    
    # SoftMask validation
    if config['softMask'] == 'yes':
        print("Softmasked genome will be generated")
    else:
        print("Softmasked genome will not be generated")
    
    # Margin validation
    if config['margin'] == 'yes':
        print("Short TE sequences <100bp will be removed from annotation")
    else:
        print("Short TE sequences <100bp will not be removed from annotation")
    
    # Helitron validation
    if config['heli'] == 'yes':
        print("HELITRON detection will be run using HELIANO")
    else:
        print("HELITRON detection will not be run")
    
    print("\nPlease cite the following paper when using this software:")
    print("Baril, T., Galbraith, J. and Hayward, A., 2024. Earl Grey: a fully automated user-friendly transposable element annotation and analysis pipeline. Molecular Biology and Evolution, 41(4), p.msae068.")
    
    return config

def check_script_directories(script_dir):
    """Check if required script directories exist"""
    if not os.path.isdir(script_dir):
        print("ERROR: Script directory variable not set, please run the configure script in the Earl Grey directory before attempting to run Earl Grey")
        sys.exit(1)
    
    testrainer_dir = os.path.join(script_dir, "TEstrainer")
    if not os.path.isdir(testrainer_dir):
        print("ERROR: teStrainer module not found, please check all modules are present and run the configure script in the Earl Grey directory before attempting to run Earl Grey")
        sys.exit(1)

def check_biocontainer():
    """Check and configure biocontainer installation if needed"""
    repeatmasker_lib = "/usr/local/share/RepeatMasker/Libraries/"
    dfam_file = os.path.join(repeatmasker_lib, "Dfam.h5")
    
    if os.path.isdir(repeatmasker_lib) and os.path.isfile(dfam_file):
        try:
            with open(dfam_file, 'r') as f:
                content = f.read()
                if "Placeholder" in content:
                    response = input("Are you using the biocontainer installation? If yes, we can configure RepeatMasker to work for you:[YyNn]")
                    if response.lower() in ['y', 'yes']:
                        configure_biocontainer(repeatmasker_lib)
                    else:
                        sys.exit(0)
        except Exception as e:
            print(f"Warning: Could not read Dfam.h5 file: {e}")

def configure_biocontainer(lib_path):
    """Configure biocontainer RepeatMasker installation"""
    try:
        os.chdir(lib_path)
        
        # Download Dfam database
        print("Downloading Dfam database...")
        urllib.request.urlretrieve(
            "https://dfam.org/releases/Dfam_3.7/families/Dfam_curatedonly.h5.gz",
            "Dfam_curatedonly.h5.gz"
        )
        
        # Extract and replace
        subprocess.run(["gunzip", "Dfam_curatedonly.h5.gz"], check=True)
        if os.path.exists("Dfam.h5"):
            os.rename("Dfam.h5", "Dfam.h5.bak")
        os.rename("Dfam_curatedonly.h5", "Dfam.h5")
        
        # Configure RepeatMasker
        os.chdir("/usr/local/share/RepeatMasker/")
        subprocess.run([
            "perl", "configure",
            "-rmblast_dir=/usr/local/bin",
            "-libdir=/usr/local/share/RepeatMasker/Libraries",
            "-trf_prgm=/usr/local/bin/trf",
            "-default_search_engine=rmblast"
        ], check=True)
        
        print("Biocontainer configuration completed successfully")
        
    except Exception as e:
        print(f"Error configuring biocontainer: {e}")
        sys.exit(1)


def check_dfam39():
    """Check for Dfam 3.9 configuration requirements"""
    try:
        # Get RepeatMasker path
        rm_path = subprocess.run(["which", "RepeatMasker"], 
                                capture_output=True, text=True, check=True).stdout.strip()
        library_path = rm_path.replace("/bin/RepeatMasker", "/share/RepeatMasker/Libraries/famdb")
        
        expected_file = os.path.join(library_path, "dfam39_full.0.h5")
        config_file = os.path.join(library_path, "rmlib.config")
        complete_marker = os.path.join(library_path, ".earlgrey.config.complete")
        
        # Check if configuration is needed
        #Check with Toby should these be "OR" instead of "AND"? Should there be a more complex check?
        if (os.path.isdir(library_path) and 
            len([f for f in os.listdir(library_path) if os.path.isfile(os.path.join(library_path, f))]) == 2 and
            os.path.isfile(expected_file) and 
            os.path.isfile(config_file) and
            not os.path.exists(complete_marker)):
            
            generate_dfam39_config_script(library_path, rm_path)
            sys.exit(2)
        #Check with Toby for the following
        #else :
        #    open(complete_marker, "w").close()
        return complete_marker
            
    except subprocess.CalledProcessError:
        print("Warning: Could not locate RepeatMasker installation")
    except Exception as e:
        print(f"Warning: Error checking Dfam 3.9 configuration: {e}")

def generate_dfam39_config_script(library_path, rm_path):
    """Generate configuration script for Dfam 3.9"""
    script_path = os.path.join(os.getcwd(), "configure_dfam39.sh")
    
    print("WARNING: Earl Grey v6.1.0 has updated to Dfam v3.9.")
    print("Before using Earl Grey, you MUST download the required partitions from Dfam")
    print(f"Configuration script generated: {script_path}")
    
    script_content = f"""#!/bin/bash
    # Configuration script for Dfam 3.9

    # Change directory to the famdb library location
    cd {library_path}/

    # Download partitions (modify the range [0-16] as needed)
    curl -o 'dfam39_full.#1.h5.gz' 'https://dfam.org/releases/current/families/FamDB/dfam39_full.[0-16].h5.gz'

    # Decompress files
    gunzip *.gz

    # Move to RepeatMasker directory
    cd {rm_path.replace('/bin/RepeatMasker', '/share/RepeatMasker/')}

    # Backup existing files
    mv {library_path}/min_init.0.h5 {library_path}/min_init.0.h5.bak

    # Reconfigure RepeatMasker
    perl ./configure \\
        -libdir {library_path.replace('/famdb', '')} \\
        -trf_prgm {rm_path.replace('/bin/RepeatMasker', '/bin/trf')} \\
        -rmblast_dir {rm_path.replace('/bin/RepeatMasker', '/bin')} \\
        -hmmer_dir {rm_path.replace('/bin/RepeatMasker', '/bin')} \\
        -abblast_dir {rm_path.replace('/bin/RepeatMasker', '/bin')} \\
        -crossmatch_dir {rm_path.replace('/bin/RepeatMasker', '/bin')} \\
        -default_search_engine rmblast

    # Mark configuration as complete
    touch {library_path}/.earlgrey.config.complete
    """
    
    with open(script_path, 'w') as f:
        f.write(script_content)
    
    #os.chmod(script_path, 0o755)
    print(f"Make the script executable and run: chmod +x {script_path} && ./{script_path}")


def make_directories(directory, species, RepSpec=None, startCust=None, heli=None):
    outdir = os.path.join(directory, f"{species}_EarlGrey")
    os.makedirs(outdir, exist_ok=True)

    if RepSpec or startCust:
        os.makedirs(os.path.join(outdir, f"{species}_RepeatMasker"), exist_ok=True)
    os.makedirs(os.path.join(outdir, f"{species}_Database"), exist_ok=True)
    os.makedirs(os.path.join(outdir, f"{species}_RepeatModeler"), exist_ok=True)
    os.makedirs(os.path.join(outdir, f"{species}_strainer"), exist_ok=True)
    os.makedirs(os.path.join(outdir, f"{species}_RepeatMasker_Against_Custom_Library"), exist_ok=True)
    os.makedirs(os.path.join(outdir, f"{species}_RepeatLandscape"), exist_ok=True)
    os.makedirs(os.path.join(outdir, f"{species}_mergedRepeats"), exist_ok=True)
    os.makedirs(os.path.join(outdir, f"{species}_summaryFiles"), exist_ok=True)
    if heli:
        os.makedirs(os.path.join(outdir, f"{species}_heliano"), exist_ok=True)
    return outdir