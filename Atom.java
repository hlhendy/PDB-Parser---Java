public class Atom{
  String name;
  Coords coordinates;
  
  //CONSTRUCTORS
  public Atom(){
    this.name = new String();
    this.coordinates = new Coords();
  }
  
  public Atom(String name){
    this.name = name;
    this.coordinates = new Coords();
  }

  public Atom(String name, Coords c){
    this.name = name;
    this.coordinates = c;
  }
  
  //METHODS///////////////////////////////            
  ////getName();
  ////getCoords();
  ////setCoords();
  ////toString();
  ////////////////////////////////////////
  public String getName(){
    return name;
  }
  
  public Coords getCoords(){
    return coordinates;
  }
  
  public void setCoords(Coords c){
    this.coordinates = c;
  }
  
  public void setCoords(double x, double y, double z){
    Coords c = new Coords(x,y,z);
    this.coordinates = c;
  }
  
  @Override
  public String toString(){
    return String.format("%s: %s\n", name, coordinates.toString());
  }
}