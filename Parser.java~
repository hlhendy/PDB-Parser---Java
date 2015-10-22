import java.util.*;
import java.io.*;

public class Parser{
  //PDB file parser
  public static Structure PDB(String filename){
    Model currModel = new Model();
    Atom currAtom = new Atom();
    Structure currStructure = new Structure(filename);
    int numModels = 0;
    try{
    Scanner sc = new Scanner(new File(filename));
    while(sc.hasNextLine()){
      String line = sc.nextLine();
      Scanner lineScanner = new Scanner(line);
      while(lineScanner.hasNext()){
        String next = lineScanner.next();
        //if token is MODEL - create new model and add to structure
        if(next.equals("MODEL") && lineScanner.hasNext()){
          //read name and create model
          if(numModels % 50 == 0){
            System.out.println("Another 50 completed\n");
          }
          String modelName = String.format("MODEL %s", numModels);
          currModel = new Model(modelName);
          numModels++;
        }
        //if token is ATOM - create new atom in current model
        if(next.equals("ATOM")){
          lineScanner.next();
          String atomName = lineScanner.next();
          //skip unwanted info
          for(int i=0; i<3; i++){
            lineScanner.next();
          }
          Double x = Double.parseDouble(lineScanner.next());
          Double y = Double.parseDouble(lineScanner.next());
          Double z = Double.parseDouble(lineScanner.next());
          currAtom = new Atom(atomName, new Coords(x, y, z));
          currModel.addAtom(currAtom);
        }
        //if token is TER - end model
//            if(next.equals("TER")){
//              currStructure.addModel(currModel);
//            }
        //if token is END - end model 
        if(next.equals("END")){
          currStructure.addModel(currModel);
        }
        //if token is ENDMDL - end model
//               if(next.equals("ENDMDL"){
//                 currStructure.addModel(currModel);
//               }
      }
      lineScanner.close();
    }
    sc.close();
    }
    catch(FileNotFoundException e){
     System.out.print("Please enter a valide file name.");
     return new Structure();
    }
    return currStructure;
  }
}