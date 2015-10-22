import java.util.*;
import java.io.*;
//CreateContactMap
public class CreateContactMap{
  public static void main(String[] args){
    Structure structure = new Structure();
    //PARSE FILE
    //check to see if filename is valid
    if(args.length != 1 || args[0] == null){
      System.out.println("Please enter a valid filename");
    }else{
      structure = Parser.PDB(args[0]);
    }
    System.out.println("Structure Created Successfully\n" + structure.toString());
    System.out.println("Typical CA count: " + structure.getModel(0).getAlphaCarbons().length);
    /////////////////CONTACT MAPS/////////////////////////
    //Binary
    //// 1: Euclidean distance between CA atoms is <= 8A
    //// 0: Distance is > 8A
    //Real-Valued
    ////Stores actual CA-CA distances
    /////////////////////////////////////////////////////
    System.out.println("Now creating contact maps\n");
    //create vectors
    //loop through vector and use distance methods in Coord class to calculate distances
    //put into real-valued map, 
    //then use criteria to determine whether to put 1 or 0 into binary map
  }
}