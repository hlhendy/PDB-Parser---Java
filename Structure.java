import java.util.*;
//Creates a structure for a list of conformations for a given protein
//A Structure has a list of models
public class Structure{
  String name;
  ArrayList<Model> modelsList;
  
  //CONSTRUCTORS
  public Structure(){
    this.name = new String();
    this.modelsList = new ArrayList<Model>();
  }
  public Structure(String name){
    this.name = name;
    this.modelsList = new ArrayList<Model>();
  }
  public Structure(String name, Model[] modelsList){
    this.name = name;
    this.modelsList = new ArrayList<Model>();
    for(int i=0; i<modelsList.length; i++){
      this.modelsList.add(modelsList[i]);
    }
  }
  
  //METHODS///////////////////////////////
  ////addModel(Model m);             
  ////getName();
  ////getModelsList();
  ////getModel(int i);
  ////getModel(String modelName);
  ////size();
  ////toString();
  ////////////////////////////////////////
  public void addModel(Model m){
    modelsList.add(m);
  }
  
  public String getName(){
    return name;
  }
  
  public ArrayList<Model> getModelsList(){
    return modelsList;
  }
  
  public Model getModel(int i){
    return modelsList.get(i);
  }
  
  public int getModel(String modelName){
    for(int i=0; i<modelsList.size(); i++){
      if(modelsList.get(i).getName().equals(modelName)){
        return i;
      }
    }
    return -1;
  }
    
  public int size(){
    return modelsList.size();
  }
  
  @Override
  public String toString(){
    StringBuilder struct = new StringBuilder();
    struct.append(String.format("%s has %d Models\n", name, modelsList.size()));
//    if(modelsList.size() > 0){
//      for(int i=0; i<modelsList.size(); i++){
//        struct.append(String.format("%s %d Atoms\n", modelsList.get(i).getName(), modelsList.get(i).size()));
//      }
//    }else{
//      struct.append("0 models");
//    }
    return struct.toString();
  }

    
}