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
  ////centroids();
  ////realign();
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
  public double lrmsd(Model m1){
    Atom[] conf1_prime = this.realign();
    Atom[] conf2_prime = m1.realign();
    
    //calculate optimal rotation
    return 0.0;
  }
  
  ////centroids(); **Just Alpha Carbons**
  public Coords centroid(){
    double[] center = {0,0,0};
    double length = this.getAlphaCarbons().length;
    for(int i=0; i < length; i++){
      center[0] += this.getAlphaCarbons()[i].getCoords().getX();
      center[1] += this.getAlphaCarbons()[i].getCoords().getY();
      center[2] += this.getAlphaCarbons()[i].getCoords().getZ();
    }
    for(int j=0; j<3; j++){
      center[j] /= length;
    }
    return new Coords(center[0], center[1], center[2]);
  }
  ////realign();
  //Get each Atom's coords, subtract centroid, update coords
  public Atom[] realign (){
    double avgX = this.centroid().getX();
    double avgY = this.centroid().getY();
    double avgZ = this.centroid().getZ();
    Atom[] caRealigned = this.getAlphaCarbons();
    for(int i=0; i<caRealigned.length; i++){
      double x = caRealigned[i].getCoords().getX();
      double y = caRealigned[i].getCoords().getY();
      double z = caRealigned[i].getCoords().getZ();
      
      caRealigned[i].setCoords((x-avgX), (y-avgY), (z-avgZ));
    }
    return caRealigned;
  }
  
  ////euclideanDistance(m1); *Array of atoms*
  public double euclideanDistance(Model m1){
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