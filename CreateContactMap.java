import java.util.*;
import java.io.*;
//CreateContactMap
public class CreateContactMap{
  
  public static void main(String[] args) throws IOException{
    int RES_DISTANCE = 4; //Atoms apart
    int BIN_CRITERIA = 8; //Angstroms
    int LRMSD_CRITERIA = 2; //Angstroms
    Structure structure = new Structure();
    Structure nativeStructure = new Structure();
    
    //PARSE FILES
    //check to see if filenames are valid
    if(args.length != 2 || args[0] == null || args[1] == null){
      System.out.println("Please enter valid filenames: <conformations>.pdb <native>.pdb");
      return;
    }else{
      nativeStructure = Parser.PDB(args[1]);
      structure = Parser.PDB(args[0]);
    }
    if(structure.getModelsList().size() >0){
      System.out.println("Structure Created Successfully\n" + structure.toString());
      System.out.println("Native Structure Created Successfully\n" + nativeStructure.toString());
    }else{
      System.out.println("Warning: PROGRAM STOPPING - NO MODELS CREATED!");
      return;
    }
    
    //////////////////////////////////////////////////////////
    //     CALCULATE lRMSD AND SORT BY LRMSD_CRITERIA       //
    //////////////////////////////////////////////////////////
    //Check to see if both models have same number of atoms -- need to output list of positions from lcs
    if(currModel.size() != nativeModel.size()){
        Model[] longestCommonSeq = lcs(currModel, nativeModel);
        currModel = longestCommonSeq[0];
        nativeModel = longestCommonSeq[1];
      }
    ArrayList<Model> withinlRMSD = new ArrayList<Model>();
    ArrayList<Model> morethanlRMSD = new ArrayList<Model>();
    ArrayList<Model> modelsList = structure.getModelsList();
    Model nativeModel = nativeStructure.getModelsList().get(0);
    for(int i=0; i <modelsList.size(); i++){
      Model currModel = modelsList.get(i);
      //then set the current model in the list to that rmsd value using currModel instead
      modelsList.get(i).setlRMSD(Distance.lrmsd(currModel, nativeModel));
      if(modelsList.get(i).getlRMSD() > LRMSD_CRITERIA){
        //sort into morethan list
        morethanlRMSD.add(modelsList.get(i));
      }else{
        //sort into within list
        withinlRMSD.add(modelsList.get(i));
      }
    }
    
    /////////////////CONTACT MAPS/////////////////////////
    //FOR EACH MODEL IN EACH SORTED LIST (AL withinlRMSD/ AL morethanlRMSD):
    //Binary
    //// 1: Euclidean distance between CA atoms is <= 8A
    //// 0: Distance is > 8A
    //Real-Valued
    ////Stores actual CA-CA distances
    /////////////////////////////////////////////////////
    System.out.println("Now creating Contact Maps.\n");
    //String output = String.format("%s_ContactMaps.txt", structure.getName());
    //PrintWriter writer = new PrintWriter(new File(output));
    ContactMap[] within_realValueMaps = new ContactMap[structure.size()];
    ContactMap[] within_binaryMaps =  new ContactMap[structure.size()];
    ContactMap[] morethan_realValueMaps = new ContactMap[structure.size()];
    ContactMap[] morethan_binaryMaps =  new ContactMap[structure.size()];
    ////////////////////////////////////////////////
    //
    //             WITHIN lRMSD CRITERIA
    //
    ////////////////////////////////////////////////
    for(int i=0; i<withinlRMSD.size(); i++){
      //create atom vectors for current model
      Atom[] ca_within = withinlRMSD.get(i).getAlphaCarbons();
      //DO NOT NEED TO BE MATRICES -- PUT IN WEKA FORMAT -- vectors
      double[][] within_realValueMap = new double[ca_within.length][ca_within.length];
      int[][] within_binMap = new int[ca_within.length][ca_within.length];
      try{
        //CALC DISTANCES FOR ALPHA CARBONS (within(k) AND ADD TO MAPS
        //for each alpha carbon
        for(int k=0; k<(ca_within.length-RES_DISTANCE); k++){
          //for each alpha carbon 4 away from current alpha carbon (within(n) and morethan(p))
          for(int n=k+RES_DISTANCE; n<ca_within.length; n++){
            within_realValueMap[k][n] = Distance.euclideanDistance(ca_within[k], ca_within[n]);
//            writer.printf(" %s  |  %s  |  %.2f  \n", 
//                                       alpha_carbons[j].toString(), alpha_carbons[k].toString(), realValueMap[j][k]);
            //add to binary maps
            if(within_realValueMap[k][n] > BIN_CRITERIA){
              within_binMap[k][n] = 0;
            }else{
              within_binMap[k][n] = 1;
            }
          }
        }
      }
      catch(Exception e){
        System.out.print("WARNING! There was an error creating the contact map.\n");
      }
      //ADD MAP TO LISTS OF MAPS
      ContactMap within_Temp = new ContactMap(String.format("Model %d", i), within_realValueMap, ca_within);
      
      //PRINT MAPS TO FILE -- CHANGE TO USE ContactMaps's toString()
//      writer.println(String.format("REAL VALUE MAP %d", (i+1)));
//      for(int m=0; m<realValueMap.length; m++){
//        for(int n=0; n<realValueMap[0].length; n++){
//          writer.print(String.format("%.2f ", realValueMap[m][n]));
//        }
//        writer.println("");
//      }
    }
  //////////////////////////////////////////////////////////////
  //
  //             MORE THAN lRMSD_CRITERIA
  //
  //////////////////////////////////////////////////////////////
  for(int j=0; j<morethanlRMSD.size(); j++){
      //create atom vectors for current model
      Atom[] ca_morethan = morethanlRMSD.get(j).getAlphaCarbons();
      //DO NOT NEED TO BE MATRICES -- PUT IN WEKA FORMAT -- vectors
      double[][] morethan_realValueMap = new double[ca_morethan.length][ca_morethan.length];
      int[][] morethan_binMap = new int[ca_morethan.length][ca_morethan.length];
      try{
        //CALC DISTANCES FOR ALPHA CARBONS morethan(m) AND ADD TO MAPS
        //for each alpha carbon
        for(int m=0; m<(ca_morethan.length-RES_DISTANCE); m++){
          //for each alpha carbon 4 away from current alpha carbon (within(n) and morethan(p))
          for(int p=m+RES_DISTANCE; p<ca_morethan.length; p++){
            morethan_realValueMap[m][p] = Distance.euclideanDistance(ca_morethan[m], ca_morethan[p]);
//            writer.printf(" %s  |  %s  |  %.2f  \n", 
//                                       alpha_carbons[j].toString(), alpha_carbons[k].toString(), realValueMap[j][k]);
            //add to binary maps
            if(morethan_realValueMap[m][p] > BIN_CRITERIA){
              morethan_binMap[m][p] = 0;
            }else{
              morethan_binMap[m][p] = 1;
            }
          }
        }
      }
      catch(Exception e){
        System.out.print("WARNING! There was an error creating the contact map.\n");
      }
      //ADD MAP TO LISTS OF MAPS
      ContactMap morethan_Temp = new ContactMap(String.format("Model %d", j), morethan_realValueMap, ca_morethan);
    }
    //writer.close();
    System.out.println("Contact Maps Created.");
  }
}