import java.util.*;
import java.io.*;
//CreateContactMap
public class CreateContactMap{
  
  public static void main(String[] args){
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
    System.out.println("Structure Created Successfully\n" + structure.toString());
    System.out.println("CA count: " + structure.getModel(0).getAlphaCarbons().length);
//    for(int i=0; i<structure.getModel(0).getAlphaCarbons().length; i++){
//      System.out.println(structure.getModel(0).getAlphaCarbons()[i]);
//    }

    /////////////////CONTACT MAPS/////////////////////////
    //Binary
    //// 1: Euclidean distance between CA atoms is <= 8A
    //// 0: Distance is > 8A
    //Real-Valued
    ////Stores actual CA-CA distances
    /////////////////////////////////////////////////////
    System.out.println("Now creating contact maps\n");
    //create vectors
    int numCA = structure.getModel(0).getAlphaCarbons().length;
    ArrayList<Model> models = structure.getModelsList();
    double[][] realValueMap = new double[numCA][numCA];
    int[][] binMap = new int[numCA][numCA];
    try{
      for(int i=0; i < (numCA-RES_DISTANCE); i++){
        for(int j=i+RES_DISTANCE; j < numCA; j++){
          realValueMap[i][j] = models.get(i).euclideanDistance(models.get(j));
          if(realValueMap[i][j] > 8){
            binMap[i][j] = 0;
          }else{
            binMap[i][j] = 1;
          }
        }
      }
    }
    catch(Exception e){
      System.out.print("WARNING! There was an error creating the contact map.\n");
    }
    
    //PRINT MAPS -- CHANGE TO WRITE TO FILE
    StringBuilder string = new StringBuilder("");
    for(int i=0; i<realValueMap.length; i++){
      for(int j=0; j<realValueMap[0].length; j++){
          string.append(String.format("%.2f ", realValueMap[i][j]));
      }
      string.append("\n");
    }
    System.out.print(string);
  }
}