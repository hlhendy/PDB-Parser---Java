import java.util.*;
import java.io.*;
//CreateContactMap
public class CreateContactMap{
  
  public static void main(String[] args) throws IOException{
    //TODO: make these part of args
    int RES_DISTANCE = 4; //Atoms apart
    int BIN_CRITERIA = 8; //Angstroms
    int LRMSD_CRITERIA = 2; //Angstroms
    Structure structure = new Structure();
    Structure nativeStructure = new Structure();
    ArrayList<Double> lrmsds = new ArrayList<Double>();
    //Atom[] atomSequence;
    
    //PARSE FILES
    //check to see if filenames are valid
    //TODO: add check to see if list is provided as second arg
    if(args.length < 2 || args.length > 3 || args[0] == null || args[1] == null){
      System.out.println("Please enter valid filenames: <conformations>.pdb <native>.pdb <optional rmsd list>.rmsd");
      return;
    }
    if(args.length >= 2){
      nativeStructure = Parser.PDB(args[1]);
      structure = Parser.PDB(args[0]);
    }
    if(args.length == 3){
      if(args[2] != null){
        lrmsds = Parser.rmsdFile(args[2]);
      }else{
        System.out.println("Please enter a valid filename for the rmsd list.");
        return;
      }
      
    }
    if((structure.getModelsList().size() >0) && (nativeStructure.getModelsList().size() > 0)){
      System.out.println("Structure Created Successfully\n" + structure.toString());
      System.out.println("Native Structure Created Successfully\n" + nativeStructure.toString());
    }else{
      System.out.println("Warning: PROGRAM STOPPING - NO MODELS CREATED!");
      return;
    }
    
    //////////////////////////////////////////////////////////
    //     CALCULATE lRMSD AND SORT BY LRMSD_CRITERIA       //
    //////////////////////////////////////////////////////////
    //If no lrmads given, check to see if both models have same number of atoms
    //need to output list of positions from lcs
    if(lrmsds == null){
      if(nativeStructure.getModel(0).size() != structure.getModel(0).size()){
        //find longest common sequence with direction
      }
      //compute lrmsds, nativeStructure is first structure
      lrmsds = Distance.lrmsd(nativeStructure, structure);
    }
    ArrayList<Model> withinlRMSD = new ArrayList<Model>();
    ArrayList<Model> morethanlRMSD = new ArrayList<Model>();
    //////////////////////////////////////////////////////////////////////////////////
    //TODO: once lrmsds calc capability updated, update this to require matching sizes
    if(lrmsds.size() != structure.size()){
      ArrayList<Double> lrmsdsShort = new ArrayList<Double>();
      for(int i=0; i<structure.size(); i++){
        lrmsdsShort.add(lrmsds.get(i));
      }
      lrmsds = lrmsdsShort;
    }
    ///////////////////////////////////////////////////////////////////////////////////
    System.out.println("Sorting by lrmsd criteria.\n");
    for(int i=0; i<lrmsds.size(); i++){
      if(lrmsds.get(i) > LRMSD_CRITERIA){
        morethanlRMSD.add(structure.getModel(i));
      }else{
        withinlRMSD.add(structure.getModel(i));
      }
    }
    System.out.printf("Lists Sorted\n Models within lrmsd criteria: %d\n Models more than lrmsd criteria: %d\n", 
                      withinlRMSD.size(), morethanlRMSD.size());
    /////////////////CONTACT MAPS/////////////////////////
    //FOR EACH MODEL IN EACH SORTED LIST (AL withinlRMSD/ AL morethanlRMSD):
    //Binary
    //// 1: Euclidean distance between CA atoms is <= 8A
    //// 0: Distance is > 8A
    //Real-Valued
    ////Stores actual CA-CA distances
    /////////////////////////////////////////////////////
    System.out.println("Now creating Contact Maps for each list.\n");
    //String output = String.format("%s_ContactMaps.txt", structure.getName());
    //PrintWriter writer = new PrintWriter(new File(output));
    int numCAs = withinlRMSD.get(0).getAlphaCarbons().length;
    double[][] within_realValueMaps = new double[withinlRMSD.size()][numCAs];
    int[][] within_binaryMaps =  new int[withinlRMSD.size()][numCAs];
    double[][] morethan_realValueMaps = new double[morethanlRMSD.size()][numCAs];
    int[][] morethan_binaryMaps =  new int[morethanlRMSD.size()][numCAs];
    //AtomPair[] morethan_atomPairs = new AtomPair[numCAs-RES_DISTANCE];
    //AtomPair[] within_atomPairs = new AtomPair[numCAs-RES_DISTANCE];
    ////////////////////////////////////////////////
    //
    //             WITHIN lRMSD CRITERIA
    //
    ////////////////////////////////////////////////
    System.out.println("Creating maps for within-lrmsd criteria list\n");
    for(int i=0; i<withinlRMSD.size(); i++){
      //create atom vectors for current model
      Atom[] ca_within = withinlRMSD.get(i).getAlphaCarbons();
      //DO NOT NEED TO BE MATRICES -- PUT IN WEKA FORMAT -- vectors
      double[] within_realValueMap = new double[ca_within.length-RES_DISTANCE];
      int[] within_binMap = new int[ca_within.length-RES_DISTANCE];
      //within_atomPairs = new AtomPair[ca_within.length-RES_DISTANCE];
      try{
        //CALC DISTANCES FOR ALPHA CARBONS (within(k) AND ADD TO MAPS
        //for each alpha carbon
        int numContacts = 0;
        for(int k=0; k<(ca_within.length-RES_DISTANCE); k++){
          //for each alpha carbon 4 away from current alpha carbon (within(n) and morethan(p))
          for(int n=k+RES_DISTANCE; n<ca_within.length; n++){
            within_realValueMap[numContacts] = Distance.euclideanDistance(ca_within[k], ca_within[n]);
            //within_atomPairs[k] = new AtomPair(ca_within[k], ca_within[n]); 
            //add to binary maps
            if(within_realValueMap[numContacts] > BIN_CRITERIA){
              within_binMap[numContacts] = 0;
            }else{
              within_binMap[numContacts] = 1;
            }
          }
          numContacts++;
        }
      }
      catch(Exception e){
        System.out.printf("WARNING! There was an error creating contact map %d.\n", i);
        return;
      }
      //ADD MAP TO LISTS OF MAPS
      within_realValueMaps[i] = within_realValueMap;
      within_binaryMaps[i] = within_binMap;
    }
//  //////////////////////////////////////////////////////////////
//  //
//  //             MORE THAN lRMSD_CRITERIA
//  //
//  //////////////////////////////////////////////////////////////
    System.out.println("Creating maps for morethan-lrmsd criteria list\n");
    for(int i=0; i<morethanlRMSD.size(); i++){
      //create atom vectors for current model
      Atom[] ca_morethan = morethanlRMSD.get(i).getAlphaCarbons();
      //DO NOT NEED TO BE MATRICES -- PUT IN WEKA FORMAT -- vectors
      double[] morethan_realValueMap = new double[ca_morethan.length-RES_DISTANCE];
      int[] morethan_binMap = new int[ca_morethan.length-RES_DISTANCE];
      try{
        //CALC DISTANCES FOR ALPHA CARBONS (within(k) AND ADD TO MAPS
        //for each alpha carbon
        int numContacts = 0;
        for(int k=0; k<(ca_morethan.length-RES_DISTANCE); k++){
          //for each alpha carbon 4 away from current alpha carbon (within(n) and morethan(p))
          for(int n=k+RES_DISTANCE; n<ca_morethan.length; n++){
            morethan_realValueMap[numContacts] = Distance.euclideanDistance(ca_morethan[k], ca_morethan[n]);
            //morethan_atomPairs[k] = new AtomPair(ca_morethan[k], ca_morethan[n]);
            //add to binary map
            if(morethan_realValueMap[numContacts] > BIN_CRITERIA){
              morethan_binMap[numContacts] = 0;
            }else{
              morethan_binMap[numContacts] = 1;
            }
          }
          numContacts++;
        }
      }
      catch(Exception e){
        System.out.printf("WARNING! There was an error creating contact map %d.\n", i);
        return;
      }
      //ADD MAP TO LISTS OF MAPS
      morethan_realValueMaps[i] = morethan_realValueMap;
      morethan_binaryMaps[i] = morethan_binMap;
    }
    System.out.println("Contact Maps Created.");
    
    //print sample maps to file for review/check
    System.out.println("Printing sample maps for models within criteria to file.");
    PrintWriter writer = new PrintWriter(new File("../Data/within_RV_sample.txt"));
    writer.append(withinlRMSD.get(0).toString() + "\r\n");
    for(int i=0; i<within_realValueMaps[0].length; i++){
     // writer.append(within_atomPairs[i] + ": " + within_realValueMaps[0][i] + "\r\n");
    }
    writer.flush();
    System.out.println("Printing sample maps for models more than criteria to file.");
    writer = new PrintWriter(new File("../Data/morethan_RV_sample.txt"));
    writer.append(morethanlRMSD.get(0).toString() + "\r\n");
    for(int i=0; i<morethan_realValueMaps[0].length; i++){
      writer.append(morethan_realValueMaps[0][i] + "\r\n");
    }
    writer.flush();
    
    //Need to put data into arff format for weka
    //http://www.cs.waikato.ac.nz/ml/weka/arff.html
    //General example: (not arff format)
    //        Atom1-Atom5 Atom5-Atom9 .... AtomN-AtomN+4 label
    // model0 distance    distance         distance      label
    //
    //so loop through each list of models and label accordingly
    //while looping, output all of the contact map info
    //for(models in withinlrmsd)
    //  print contact map info, then either 0 (morethan) or 1 (within)
    System.out.println("Creating arff files for Weka.");
    String filename = "ContactMapFeatures.arff";
    String directory = "../Data/";
    File file = new File(directory, filename);
    writer = new PrintWriter(file);
    //Comments
    writer.append("% 1. Title: Contact Map Features\n%\n%\n");
    //Header: Relation and Attributes
    writer.append("@RELATION protein-conformations\n");
    writer.append("\n");
    //Attributes: each atom-atom pair
    //for each atom pair:
    ////@ATTRIBUTE atomname-atomname NUMERIC
    ////...
    ////@ATTRIBUTE label NUMERIC
    for(int i=0; i<within_realValueMaps[0].length; i++){
      writer.append(String.format("@ATTRIBUTE pair%d NUMERIC\n", i));
    }
    writer.append("@ATTRIBUTE class NUMERIC\n");
    writer.append("\n");
    //Data
    //@DATA
    //for each conformation
    ////5.1,3.5,...,10.5,0 -------->each distance from contact map, then 0 or 1 for label
    
    writer.append("@DATA\n");
    for(int i=0; i<within_realValueMaps.length; i++){
      StringBuilder sb = new StringBuilder();
      //each value in current map separated by a comma
      for(int j=0; j<within_realValueMaps[i].length; j++){
        sb.append(within_realValueMaps[i][j]);
        sb.append(",");
      }
      //append label at end of current contact map data
      sb.append("1\n");
      writer.append(sb.toString());
    }
    for(int i=0; i<morethan_realValueMaps.length; i++){
      StringBuilder sb = new StringBuilder();
      //each value in current map separated by a comma
      for(int j=0; j<morethan_realValueMaps[i].length; j++){
        sb.append(morethan_realValueMaps[i][j]);
        sb.append(",");
      }
      //append label at end of current contact map data
      sb.append("0\n");
      writer.append(sb.toString());
    }
    System.out.printf("File Created in directory %s\n", directory);
    writer.close();
  }
}
