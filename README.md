# CRISPR Design for Parasitic Nematodes

Welcome to the CRISPR Design for Parasitic Nematodes repository. This tool facilitates the design of CRISPR experiments targeting parasitic nematodes, providing features such as off-target search, results visualization, and primer design.

## Latest Version Updates
- **1.0.9:** Added graphic of sgRNA targets.
- **1.0.8:** Added preloader.
- **1.0.7:** Integrated Primer3 output into results table.
- **1.0.6:** CRISPRATER score column added to results.
- **1.0.5:** PAMPN score column added to results.
- **1.0.4:** Results column added with PAMPN GC%.
- **1.0.3:** Added file-lock structure to balance job queue.
- **1.0.2:** New method now available to search for off-targets.
- **1.0.1:** New interface for design parameters.
- **1.0.0:** Published at Molecular & Biochemical Parasitology.

## File Structure

### Source Code
- **app.pl:** Main application script.
- **index.html:** HTML file for the main application interface.
- **run_crispr.php:** PHP script to execute the CRISPR design.

### Styles and Assets
- **/css:** Stylesheets for the application.
- **/favicon:** Icons for different platforms.
- **/fonts:** Font files used in the application.
- **/html:** HTML files related to the application's CGI scripts.
- **/img:** Images used in the application.
- **/index:** Index files related to the Bowtie tool.
- **/js:** JavaScript files used in the application.
- **/sass:** Sass files and styles related to Bootstrap.

### CGI Scripts
- **/html/cgi-bin:** CGI scripts for Primer3 and related configurations.
- **/html/cgi-bin/cache:** Cached files used by CGI scripts.
- **/html/cgi-bin/primer3_config:** Primer3 configurations and interpretations.

## Usage

1. Execute the `run_crispr.php` script to initiate the CRISPR design.
2. Access the main application interface via `index.html`.
3. Explore the results and design parameters for CRISPR experiments.

## Contributions

Contributions to the development and improvement of this CRISPR design tool are welcome. If you encounter any issues or have suggestions, please create an issue in the repository.

Happy CRISPR designing!
