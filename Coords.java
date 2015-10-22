public class Coords{
  Double x;
  Double y;
  Double z;
  
  //CONSTRUCTORS
  public Coords(){
    this.x = 0.0;
    this.y = 0.0;
    this.z = 0.0;
  }
  public Coords(Double x, Double y, Double z){
    this.x = x;
    this.y = y;
    this.z = z;
  }
  public Coords(Coords c){
    this.x = c.getX();
    this.y = c.getY();
    this.z = c.getZ();
  }
  
  //METHODS///////////////////////////////
  ////getX, getY, getZ
  ////lrmsd(Coords c1);
  ////euclideanDistance(Coords c1);
  ////manhattanDistance(Coords c1);
  ////distanceSquared(Coords c1);
  ////toString();
  ////////////////////////////////////////
  public Double getX(){
    return x;
  }
  
  public Double getY(){
    return y;
  }
  
  public Double getZ(){
    return z;
  }
  
  public double euclideanDistance(Coords c1){
    return Math.sqrt(this.distanceSquared(c1));
  }
  public double manhattanDistance(Coords c1){
    double deltaX = this.x - c1.getX();
    double deltaY = this.y - c1.getY();
    double deltaZ = this.z - c1.getZ();
    double result = Math.abs(deltaX) + Math.abs(deltaY) + Math.abs(deltaZ);
    return result;
  }
  public double distanceSquared(Coords c1){
    double deltaX = this.x - c1.getX();
    double deltaY = this.y - c1.getY();
    double deltaZ = this.z - c1.getZ();
    double result = (deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ);
    return result;
  }
  
  @Override
  public String toString(){
    return String.format("[%f, %f, %f]\n", x, y, z);
  }
}