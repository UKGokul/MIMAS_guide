#!/usr/bin/env python3
"""
Build script for MIMAS documentation website.
Parses flow.md and mimas_A_2014.f to generate JSON data files.
"""

import re
import json
import os
from pathlib import Path

# Paths
BASE_DIR = Path(__file__).parent.parent
FLOW_FILE = BASE_DIR / "flow.md"
FORTRAN_FILE = BASE_DIR / "mimas_A_2014 (1).f"
OUTPUT_DIR = BASE_DIR / "data"

def parse_flow_md():
    """Parse flow.md to extract stage information."""
    with open(FLOW_FILE, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Split by headings # <number>) or # <number>
    pattern = r'#\s*(\d+)\)[ \t]*(.*?)(?=\n#\s*\d+[ \t]*\)|\Z)'
    # Alternative pattern for headings without )
    pattern2 = r'#\s*(\d+)[ \t]+(.*?)(?=\n#\s*\d+[ \t]+|\Z)'
    
    stages = []
    
    # Split content to find sections
    lines = content.split('\n')
    current_stage = None
    current_content = []
    
    stage_pattern = re.compile(r'^#\s*(\d+)\)[ \t]*(.*)$')
    stage_pattern2 = re.compile(r'^#\s*(\d+)[ \t]+(.*)$')
    
    for line in lines:
        # Check if this is a new stage heading
        match = stage_pattern.match(line) or stage_pattern2.match(line)
        if match:
            # Save previous stage
            if current_stage is not None:
                # Extract overview from content (everything after title, excluding block_diagram lines)
                overview_lines = []
                for cl in current_content:
                    cl_stripped = cl.strip()
                    if cl_stripped and not cl_stripped.startswith('block_diagram'):
                        overview_lines.append(cl_stripped)
                
                stages.append({
                    "id": int(current_stage),
                    "title": current_title,
                    "overview": '\n'.join(overview_lines)
                })
            
            # Start new stage
            current_stage = match.group(1)
            current_title = match.group(2).strip()
            current_content = []
        else:
            if current_stage is not None:
                current_content.append(line)
    
    # Don't forget the last stage
    if current_stage is not None:
        overview_lines = []
        for cl in current_content:
            cl_stripped = cl.strip()
            if cl_stripped and not cl_stripped.startswith('block_diagram'):
                overview_lines.append(cl_stripped)
        
        stages.append({
            "id": int(current_stage),
            "title": current_title,
            "overview": '\n'.join(overview_lines)
        })
    
    return stages

def parse_fortran():
    """Parse Fortran file to extract subroutines and call graph."""
    with open(FORTRAN_FILE, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()
    
    # Find all subroutine definitions
    sub_pattern = re.compile(r'^\s*SUBROUTINE\s+(\w+)', re.IGNORECASE | re.MULTILINE)
    sub_matches = list(sub_pattern.finditer(content))
    
    subroutines = {}
    
    for match in sub_matches:
        sub_name = match.group(1).lower()
        start_pos = match.start()
        
        # Find end of subroutine (next SUBROUTINE or FUNCTION or end of file)
        end_pos = len(content)
        for other_match in sub_matches:
            if other_match.start() > start_pos:
                end_pos = other_match.start()
                break
        
        # Extract subroutine body
        sub_body = content[start_pos:end_pos]
        
        # Find all calls within this subroutine
        call_pattern = re.compile(r'\bCALL\s+(\w+)', re.IGNORECASE)
        calls = [call.group(1).lower() for call in call_pattern.finditer(sub_body)]
        
        # Get snippet (first ~30 lines or 500 chars)
        snippet_lines = sub_body.split('\n')[:30]
        snippet = '\n'.join(snippet_lines)
        if len(sub_body.split('\n')) > 30:
            snippet += '\n...'
        
        subroutines[sub_name] = {
            "name": sub_name,
            "calls": list(set(calls)),
            "snippet": snippet.strip()
        }
    
    return subroutines

def generate_science_view(stage_id, stage_overview):
    """Generate science view from flow.md content."""
    # Science view is derived from the overview - we don't invent new physics
    science_lines = []
    
    # Extract key physics concepts mentioned
    if 'advection' in stage_overview.lower() or 'transport' in stage_overview.lower():
        science_lines.append("Transport processes move water vapor and particles through the atmosphere.")
    
    if 'diffusion' in stage_overview.lower() or 'turbkzz' in stage_overview.lower():
        science_lines.append("Vertical turbulent diffusion redistributes water vapor based on atmospheric stability.")
    
    if 'photolysis' in stage_overview.lower() or 'lyman' in stage_overview.lower():
        science_lines.append("Solar UV radiation dissociates water vapor molecules, affecting H₂O concentrations.")
    
    if 'ice' in stage_overview.lower() or 'particle' in stage_overview.lower():
        science_lines.append("Ice nucleation occurs on dust particles when supersaturation conditions are met.")
    
    if 'saturation' in stage_overview.lower() or 'kelvin' in stage_overview.lower():
        science_lines.append("Saturation vapor pressure depends on temperature and particle curvature (Kelvin effect).")
    
    if 'sedimentation' in stage_overview.lower() or 'gravity' in stage_overview.lower():
        science_lines.append("Ice particles fall due to gravitational settling, removing water from upper atmosphere.")
    
    if not science_lines:
        science_lines.append("See overview for details on this model's physics processes.")
    
    return " ".join(science_lines)

def find_relevant_subroutines(stage, subroutines):
    """Find subroutines mentioned in stage text."""
    relevant = []
    stage_text = (stage['title'] + ' ' + stage['overview']).lower()
    
    for sub_name in subroutines:
        # Check if subroutine name appears in stage text
        if sub_name in stage_text:
            relevant.append(sub_name)
        # Also check for variations (sub_ prefix)
        sub_var = sub_name.replace('sub_', '')
        if sub_var in stage_text and sub_var not in ['init', 'lesedyn', 'zenit', 'photolyse', 'tabelle', 'diffu', 'transwalcek', 'tracertransp']:
            if sub_name not in relevant:
                relevant.append(sub_name)
    
    return relevant

def main():
    """Main build function."""
    print("Building MIMAS documentation data files...")
    
    # Create output directory
    OUTPUT_DIR.mkdir(exist_ok=True)
    
    # Parse flow.md
    print("Parsing flow.md...")
    stages = parse_flow_md()
    print(f"  Found {len(stages)} stages")
    
    # Parse Fortran file
    print("Parsing mimas_A_2014.f...")
    subroutines = parse_fortran()
    print(f"  Found {len(subroutines)} subroutines")
    
    # Build flow.json with enhanced data
    flow_data = []
    for stage in stages:
        relevant_subs = find_relevant_subroutines(stage, subroutines)
        science_view = generate_science_view(stage['id'], stage['overview'])
        
        # Extract mentioned variables
        variables = []
        stage_text = (stage['title'] + ' ' + stage['overview']).lower()
        common_vars = ['hm', 'hminit3d', 'um1', 'um2', 'vm1', 'vm2', 'wm1', 'wm2', 
                       'xfeld', 'yfeld', 'zfeld', 'dttrans', 'turbkzz', 'rfeld', 'rinit']
        for var in common_vars:
            if var in stage_text:
                variables.append(var)
        
        flow_data.append({
            "id": stage['id'],
            "title": stage['title'],
            "overview": stage['overview'],
            "science_view": science_view,
            "relevant_subroutines": relevant_subs if relevant_subs else ["Not identified yet"],
            "variables": variables if variables else ["Not identified yet"]
        })
    
    # Write flow.json
    flow_json_path = OUTPUT_DIR / "flow.json"
    with open(flow_json_path, 'w', encoding='utf-8') as f:
        json.dump(flow_data, f, indent=2)
    print(f"  Written: {flow_json_path}")
    
    # Write subroutines.json
    subroutines_json_path = OUTPUT_DIR / "subroutines.json"
    with open(subroutines_json_path, 'w', encoding='utf-8') as f:
        json.dump(subroutines, f, indent=2)
    print(f"  Written: {subroutines_json_path}")
    
    print("\nBuild complete!")

if __name__ == "__main__":
    main()
