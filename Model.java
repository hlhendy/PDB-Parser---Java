import java.util.ArrayList;

public class Model{
  String name;
  ArrayList<Atom> atomsList;
  
  //CONSTRUCTORS
  public Model(){
    this.name = new String();
    this.atomsList = new ArrayList<Atom>();
  }
  public Model(String name){
    this.name = name;
    this.atomsList = new ArrayList<Atom>();
  }
  public Model(String name, Atom[] atomsList){
    this.name = name;
    this.atomsList = new ArrayList<Atom>();
    for(int i=0; i<atomsList.length; i++){
      this.atomsList.add(atomsList[i]);
    }
  }
  
  //METHODS///////////////////////////////
  ////addAtom(Atom a);             
  ////getName();
  ////getAtomsList();
  ////getAlphaCarbons();
  ////size();
  ////lrmsd(Model m1);
  ////centroids(Model m1);
  ////realign(Model m1);
  ////euclideanDistance(m1); *Array of atoms*
  ////toString();
  ////////////////////////////////////////
  public void addAtom(Atom a){
    atomsList.add(a);
  }
  
  public String getName(){
    return name;
  }
  
  public ArrayList<Atom> getAtomsList(){
    return atomsList;
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
  
  ////lrmsd(Model m1);
  ////centroids(Model m1);
  ////realign(Model m1);
  ////euclideanDistance(m1); *Array of atoms*
  double euclideanDistance(Model m1){
    double distance = 0.0;
    Atom[] mCA = this.getAlphaCarbons();
    Atom[] m1CA = m1.getAlphaCarbons();
    for(int i=0; i < this.getAlphaCarbons().length; i++){
      distance = distance + mCA[i].getCoords().euclideanDistance(m1CA[i].getCoords());
    }
    return distance;
  }
  
  @Override
  public String toString(){
    return String.format("%s\n--------\n%d Atoms\n%d Alpha Carbons\n", name, size(), getAlphaCarbons().length);
  }
}