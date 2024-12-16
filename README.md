# PhyloTreeGUI

A tool for visualizing interactive phylogenetic trees. PhyloTreeGUI takes a phyloXML file as input and produces the following outputs:

1) a BED file containing the annotation of the selected HORs (Higher-Order Repeats) in the centromere sequence;
2) a file describing the composition of the selected HORs;
3) a set of files illustrating the monomers belonging to the families that form the selected HORs.

PhyloTreeGUI can be launched from a terminal:

python main.py [--path <path/to/dir/with/phylo/xml/files>] [--folder <dir/where/to/save/output/files>] [--fullscreen]

Arguments within squared brackets are optional and can be used to specify:

1) path to the directory containing the phyloXML files (default is the current directory);
2) folder/directory where output files will be stored (default is the current directory);
3) whether PhyloTreeGUI has to be run in full-screen mode (default is false).

When the GUI is opened, the user can see a drop-down menu on the upper-right corner, which shows the list of phyloXML files in the input folder that can be selected. The GUI allows to perform different actions through the following buttons:

1) LoadFile
2) PlotTree
3) ShowData
4) ZoomIn
5) ZoomOut
6) GetOutput
7) Reset

After selecting a phyloXML file from the drop-down menu, the user shall click on the “LoadFile” button in order to load all the information contained in the file (i.e., the monomer and the HOR trees). Then, the user can see the HOR tree by pressing the “PlotTree” button: it consists of several branches originating from the root and containing circles/nodes, which indicate the HORs found. The HOR tree is interactive: the nodes are clickable and can be used to select one or more HORs, which will assume a different color. Once the user has selected a set of HORs to analyze, he/she can visualize the structure of the HORs (i.e., the families of monomers composing each of them), the locations of the HORs in the chromosome sequence and the monomer tree by clicking the “ShowData” button. In particular, a new panel appears on the bottom right corner of the PhyloTreeGUI, displaying each selected HOR (identified by the same color as the clicked node) along with a schematic representation of its composition. This schematic is composed of N colored blocks, where N corresponds to the number of clades composing the HOR, and each color is specific to a clade. The same clade-specific colors also appear in the monomer tree (displayed on the left side of the PhyloTreeGUI) to distinguish all the branches involved in the HORs. An additional section in the same panel shows a grey bar representing the input sequence, with colored portions indicating the exact location of the selected HOR.
Moreover, the monomer tree might include nodes, which represent collapsed sub-trees containing more than 10 clades. In order to inspect nodes, the user can click one or more of them and then press the “ZoomIn” button: the clicked node(s) is(are) expanded and a light-gray rectangle identifies the corresponding region(s) of the tree. The resulting plot shows only part of the tree, the user can see the rest by using a scroll bar on the right side of the plot. To exit the visualization, the user shall press the “ZoomOut” button, which shows the previous monomer tree with no expanded nodes. The user can make different choices about nodes to inspect and perform a different analysis. The “Reset” button allows to reset the GUI: it clears the selection of the HORs in the HOR tree and removes all the other objects that have been visualized (e.g., the families of monomers constituting the selected HORs, their locations in the chromosome sequence and the monomer tree). The user can make a new selection of different HORs and perform a different analysis by repeating the steps previously described.
In order to save data related to the selected HORs, the user should click the “GetOutput” button that is responsible for the generation of the output files described above. 
If the user decides to load a different phyloXML file (drop-down menu), the window containing the previous HOR tree will be cleaned. Finally, the tool shows pop-up error messages when the user performs one of the following invalid operations:

1) the user clicks on the “LoadFile” button before having selected a phyloXML from the drop-down menu;
2) the user presses the “PlotTree” button before having loaded the phyloXML file through the “LoadFile” button;
3) the user clicks on the “ShowData” button before having selected one or more HORs from the HOR tree;
4) the user presses the “ZoomIn” button before selecting the nodes to inspect;
5) the user clicks on the “ZoomOut” button, but no nodes have been inspected through the “ZoomIn” button.


