================================================================================
MIMAS Model Documentation Website
================================================================================

This is a static documentation website that explains the MIMAS (Mesospheric 
Ice Model with Adaptive Sublimation) model flow step-by-step.


================================================================================
FILE STRUCTURE
================================================================================

website/
  ├── index.html              # Landing page with link to guide
  ├── pages/
  │   └── guide.html          # Main documentation page
  ├── assets/
  │   ├── css/
  │   │   └── style.css       # All styles
  │   └── js/
  │       └── app.js          # Main application logic
  ├── data/
  │   ├── flow.json           # Stage data (21 stages parsed from flow.md)
  │   └── subroutines.json   # Subroutine data from mimas_A_2014.f
  ├── scripts/
  │   └── build_data.py       # Build script (for regenerating data files)
  ├── flow.md                 # Source documentation
  ├── mimas_A_2014 (1).f      # Fortran source code
  └── README.txt             # This file


================================================================================
CONTENT
================================================================================

The website documents the 21 stages of the MIMAS model:

1. Main program starts + defines simulation window
2. Calendar setup
3. Global initialization
4. Read initial 1D H2O profile
5. Read dynamics from LIMA
6. Numerical safety check (CFL)
7. Initialize 40 million dust particles
8. Restart/output bootstrapping
9. Precompute microphysics lookup tables
10. Loop over season → days → hours
11. Map initial H2O onto 3D grid
12. Compute solar geometry and photolysis
13. Update background dynamics
14. Eulerian H2O transport (advection)
15. Vertical diffusion of H2O
16. Photolysis again
17. Particle microphysics + transport
18. Dust reallocation
19. Write diagnostics
20. Write particle positions
21. Advance time


================================================================================
WEBSITE FEATURES
================================================================================

- Landing page with quick access to model guide
- 21-stage flowchart on the left sidebar (21 numbered stages)
- Three accordion panels for each stage:
  A) Overview - what happens in this stage (parsed from flow.md)
  B) Code View - relevant subroutines and variables
  C) Science View - physics explanation
- Clickable subroutine names showing Fortran code snippets
- Search box to filter stages by title or subroutine name
- Keyboard navigation (arrow keys)
- Mobile-responsive design (dropdown selector on mobile)
- Modal popup for viewing subroutine definitions


================================================================================
UPDATING CONTENT
================================================================================

To update the model flow documentation:

1. Edit flow.md with your changes
   - Stages are marked with "# <number>)" or "# <number>" headings
   - Add new sections as needed

2. Regenerate data files:
   python scripts/build_data.py
   
   OR manually update data/flow.json if Python is not available

3. Refresh your browser to see changes


================================================================================
RUNNING THE BUILD SCRIPT (OPTIONAL)
================================================================================

The build script parses flow.md and mimas_A_2014.f to regenerate data files:

cd c:/Users/gu01/Desktop/website
python scripts/build_data.py

Requirements:
- Python 3.6+
- Access to flow.md and mimas_A_2014.f

Output:
- data/flow.json (21 stages with overview, science view, subroutines, variables)
- data/subroutines.json (subroutine definitions, call graph, code snippets)


================================================================================
PREVIEW THE WEBSITE
================================================================================

Option 1: Using Python's http.server
--------------------------------------
cd c:/Users/gu01/Desktop/website
python -m http.server 8000

Then open: http://localhost:8000

Option 2: Using VS Code Live Server
--------------------------------------
If using VS Code:
1. Install "Live Server" extension
2. Right-click index.html → "Open with Live Server"

Option 3: Direct file access
--------------------------------------
Open index.html directly in a browser (some features may not work due to CORS)


================================================================================
TROUBLESHOOTING
================================================================================

Q: The page shows "Failed to load data files"
A: Check that data/flow.json and data/subroutines.json exist.
   If not, run: python scripts/build_data.py
   Or manually copy the JSON structure.

Q: Changes to flow.md don't appear
A: Regenerate data files and refresh your browser.

Q: Subroutines show "Not identified yet"
A: The subroutine name must appear in the stage text in flow.md.
   Or manually add the mapping in data/flow.json.

Q: The server shows directory listing instead of the page
A: Make sure you're accessing the correct URL.
   Use http://localhost:8000/ for Python server.


================================================================================
BROWSER COMPATIBILITY
================================================================================

Tested on:
- Chrome 90+
- Firefox 88+
- Safari 14+
- Edge 90+

Requires modern JavaScript features (ES6).


================================================================================
DATA FILES FORMAT
================================================================================

flow.json:
[
  {
    "id": 1,
    "title": "Stage title",
    "overview": "Full description from flow.md",
    "science_view": "Physics explanation",
    "relevant_subroutines": ["sub_name1", "sub_name2"],
    "variables": ["var1", "var2"]
  },
  ...
]

subroutines.json:
{
  "sub_name": {
    "name": "sub_name",
    "calls": ["called_sub1", "called_sub2"],
    "snippet": "First ~30 lines of subroutine"
  },
  ...
}


================================================================================
