4/6/2017


QUICK START INSTRUCTIONS
========================

These instructions follow the typical expected usage of the pair-guides script. It is assumed that the latest version of the pair-guides script has been configured. These instructions will generate candidate guide pairs for grk5, but they apply the same way to any gene.

1. Create a new gene folder. Title it "grk5"

2. Search grk5 in ChopChop (http://chopchop.cbu.uib.no/). From the graphical results page, find the drop-down menu titled "Download results" and choose "Results table". This will load the ChopChop results as a text file in a new tab.

3. Save the ChopChop results. Right-click on the text results page and select "Download" if possible. Otherwise search for a "Save as" option. Save the results to the gene folder.

4. Create a settings file. The simplest way to do this is to copy the settings file from the project directory's demo folder (likely Desktop/Scripts/pair-guides/demo/gene_block_settings.inp). The other option is to simply create a new text file named "gene_block_settings.inp" in the gene folder.

5. Edit the settings file with TextEdit. Any line that starts with a hash symbol (#) is treated as a comment and ignored. There must be three lines at minimum, formatted like this:

  CCDS_ID     <insert the CCDS id of the gene>
  input_file  <insert the name of the ChopChop results.txt file>
  strand      <insert the strand with the coding sequence ('+' or '-')>

Each line contains a key and a value. Find a sample settings file in the "demo" folder of the project directory (likely Desktop/Scripts/pair-guides/demo). This demo directory contains inputs for the grk5 gene.

6. Open the gene folder in Terminal. This can be done by dragging the gene folder to the Terminal icon in the tray or by right-clicking on the gene folder icon and selecting Services > "New Terminal at Folder".

7. Run pair-guides. This is done by typing "pair-guides" into the Terminal and hitting the Enter/Return key.

8. After the script finished running, the gene folder will contain an additional file that ends with "_pairs.csv". This is a spreadsheet containing the candidate guide pairs, some data to help evaluate these pairs, and full gene block sequences for each pair.




Email Nate Herrmann (naherrmann@gmail.com) with notes on any bugs or crashes. Include as much detail as possible.