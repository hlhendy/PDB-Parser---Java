import java.util.*;
import java.io.*;
//CreateContactMap
public class CreateContactMap{
  
  public static void main(String[] args) throws IOException{
    int RES_DISTANCE = 4;
    int BIN_CRITERIA = 8; //A
    Structure structure = new Structure();
    //PARSE FILE
    //check to see if filename is valid
    if(args.length != 1 || args[0] == null){
      System.out.println("Please enter a valid filename");
    }else{
      structure = Parser.PDB(args[0]);
    }
    String output = String.format("%s_ContactMaps.txt", structure.getName());
    System.out.println(output);
    PrintWriter writer = new PrintWriter(new File(output));
    if(structure.getModelsList().size() >0){
      System.out.println("Structure Created Successfully\n" + structure.toString());
      writer.println("CA count: " + structure.getModel(0).getAlphaCarbons().length);
//      for(int i=0; i<structure.getModel(0).getAlphaCarbons().length; i++){
//        System.out.println(structure.getModel(0).getAlphaCarbons()[i]);
//      }
    }else{
      System.out.println("Warning: PROGRAM STOPPING - NO MODELS CREATED!");
      return;
    }
    
    /////////////////CONTACT MAPS/////////////////////////
    //FOR EACH MODEL:
    //Binary
    //// 1: Euclidean distance between CA atoms is <= 8A
    //// 0: Distance is > 8A
    //Real-Valued
    ////Stores actual CA-CA distances
    /////////////////////////////////////////////////////
    System.out.println("Now creating contact maps\n");
    ArrayList<Model> models = structure.getModelsList();
    //MODELS (i)
    for(int i=0; i<models.size(); i++){
      //create vectors
      Atom[] alpha_carbons = models.get(i).getAlphaCarbons();
      double[][] realValueMap = new double[alpha_carbons.length][alpha_carbons.length];
      int[][] binMap = new int[alpha_carbons.length][alpha_carbons.length];
      try{
        //ALPHA CARBONS (j and k)
        for(int j=0; j < (alpha_carbons.length-RES_DISTANCE); j++){
          for(int k=j+RES_DISTANCE; k < alpha_carbons.length; k++){
            realValueMap[j][k] = Distance.euclideanDistance(alpha_carbons[j], alpha_carbons[k]);
            if(realValueMap[j][k] > 8){
              binMap[j][k] = 0;
            }else{
              binMap[j][k] = 1;
            }
          }
        }
      }
      catch(Exception e){
        System.out.print("WARNING! There was an error creating the contact map.\n");
      }
      //PRINT MAPS TO FILE
      writer.println(String.format("REAL VALUE MAP %d", (i+1)));
      for(int m=0; m<realValueMap.length; m++){
        for(int n=0; n<realValueMap[0].length; n++){
          writer.print(String.format("%.2f ", realValueMap[m][n]));
        }
        writer.println("");
      }
    }
    writer.close();
    System.out.println("Contact Maps Created.");
  }
}