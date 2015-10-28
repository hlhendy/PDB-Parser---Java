public class Atom{
  String name;
  double[] coords;
  
  //CONSTRUCTORS
  public Atom(){
    this.name = new String();
    coords = new double[3];
  }
  
  public Atom(String name){
    this();
    this.name = name;
  }

  public Atom(String name, double x, double y, double z){
    this(name);
    coords[0] = x;
    coords[1] = y;
    coords[2] = z;
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
  
  public double[] getCoords(){
    return coords;
  }
  
  public void setCoords(double x, double y, double z){
    coords[0] = x;
    coords[1] = y;
    coords[2] = z;
  }
  
  @Override
  public String toString(){
    StringBuilder c = new StringBuilder("");
    c.append(String.format("[%.2f, %.2f, %.2f]", coords[0], coords[1], coords[2]));
    return(String.format("%s : %s\n", name, c.toString()));
  }
}