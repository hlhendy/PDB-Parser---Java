import java.util.*;
import java.io.*;

public class Parser{
  //PDB file parser
  public static Structure PDB(String filename){
    Model currModel = new Model();
    Atom currAtom = new Atom();
    boolean inModel = false;
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
        if((next.equals("MODEL") || (next.equals("ATOM") && !inModel) && lineScanner.hasNext())){
          //read name and create model
          //int modelNum = Integer.parseInt(line.substring(5,6));
          String modelName = String.format("MODEL %s", numModels);
          currModel = new Model(modelName);
          numModels++;
          inModel = true;
        }
        //if token is ATOM and a model has been started - create new atom in current model
        if(next.equals("ATOM") && inModel){
          lineScanner.next();
          String atomName = lineScanner.next();
          //skip unwanted info
          for(int i=0; i<3; i++){
            lineScanner.next();
          }
          Double x = Double.parseDouble(lineScanner.next());
          Double y = Double.parseDouble(lineScanner.next());
          Double z = Double.parseDouble(lineScanner.next());
          currAtom = new Atom(atomName, x, y, z);
          currModel.addAtom(currAtom);
        }
        //if token is TER - end model
//            if(next.equals("TER")){
//              currStructure.addModel(currModel);
//            }
        //if token is END - end model 
        if(next.equals("END")){
          currStructure.addModel(currModel);
          inModel = false;
          if(currStructure.size() % 100 == 0){
            System.out.printf("%d models completed\n", currStructure.size());
          }
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
     System.out.print("Please enter a valid file name.");
     return new Structure();
    }
    return currStructure;
  }
  //rmsd file parser -- must  be list of doubles
  public static ArrayList<Double> rmsdFile(String filename){
    ArrayList<Double> rmsd = new ArrayList<Double>();
    try{
      Scanner sc = new Scanner(new File(filename));
      double next = 0.0;
      while(sc.hasNextDouble()){
        next = sc.nextDouble();
        rmsd.add(next);
      } 
    }
    catch(FileNotFoundException e){
      System.out.print("Please enter a valid file name.");
      return new ArrayList<Double>();
    }
    return rmsd;
  }
  
}