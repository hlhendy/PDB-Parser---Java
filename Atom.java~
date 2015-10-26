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
  ////toString();
  ////////////////////////////////////////
  public String getName(){
    return name;
  }
  
  public Coords getCoords(){
    return coordinates;
  }
  
  @Override
  public String toString(){
    return String.format("%s: %s\n", name, coordinates.toString());
  }
}