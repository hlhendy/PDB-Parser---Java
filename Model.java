import java.util.ArrayList;

public class Model{
  String name;
  double lRMSD;
  ArrayList<Atom> atomsList;
  
  //CONSTRUCTORS
  public Model(){
    this.name = new String();
    this.lRMSD = 0.0;
    this.atomsList = new ArrayList<Atom>();
  }
  public Model(String name){
    this.name = name;
    this.lRMSD = 0.0;
    this.atomsList = new ArrayList<Atom>();
  }
  public Model(String name, double lRMSD, Atom[] atomsList){
    this.name = name;
    this.lRMSD = lRMSD;
    this.atomsList = new ArrayList<Atom>();
    for(int i=0; i<atomsList.length; i++){
      this.atomsList.add(atomsList[i]);
    }
  }
  
  //METHODS///////////////////////////////
  ////addAtom(Atom a);             
  ////getName();
  ////getAtomsList();
  //getAtomNames();
  ////getAlphaCarbons();
  ////size();
  ////toString();
  ////////////////////////////////////////
  public void addAtom(Atom a){
    atomsList.add(a);
  }
  
  public String getName(){
    return name;
  }
  
  public void setlRMSD(double lRMSD){
    this.lRMSD = lRMSD;
  }
  
  public double getlRMSD(){
    return lRMSD;
  }
  
  public ArrayList<Atom> getAtomsList(){
    return atomsList;
  }
  
  public char[] getAtomNames(){
    char[] atomNames = new char[atomsList.size()];
    for(int i=0; i<atomsList.size(); i++){
      atomNames = atomsList.get(i).getName();
    }
  }
  
  public Atom[] getAlphaCarbons(){
    ArrayList<Atom> alphaCarbons = new ArrayList<Atom>();
    for(int i=0; i<atomsList.size(); i++){
      if(atomsList.get(i).getName().equals("CA")){
        alphaCarbons.add(atomsList.get(i));
      }
    }
    return alphaCarbons.toArray(new Atom[alphaCarbons.size()]);
  }
    
  public int size(){
    return atomsList.size();
  }
  
  @Override
  public String toString(){
    return String.format("%s\n--------\n%d Atoms\n%d Alpha Carbons\n", name, size(), getAlphaCarbons().length);
  }
}