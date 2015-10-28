public class ContactMap<T>{
  double[][] rvMap;
  //int[] binMap;
  String name;
  Atom[] alphaCarbons;
  
  public <T> ContactMap(String name, double[][] map, Atom[] ac){
    //Name should include protein name, model number, and type of map
    this.name = name;
    this.rvMap = map;
    this.alphaCarbons = ac;
  }
  
  public double[][] getMap(){
    return rvMap;
  }
  
  public void setMap(double[][] map){
    this.rvMap = map;
  }
  
  public String getName(){
    return name;
  }
  
  public void setName(String name){
    this.name = name;
  }
  
  public Atom[] getAlphaCarbons(){
    return alphaCarbons;
  }
  
  public void setAlphaCarbons(Atom[] ac){
    this.alphaCarbons = ac;
  }
  
  @Override
  public String toString(){
    StringBuilder sb = new StringBuilder(String.format("%s\n", name));
      for(int i=0; i<rvMap.length; i++){
        for(int j=0; i<rvMap[0].length; j++){
          sb.append(String.format("%.2f ", rvMap[i][j]));
        }
        sb.append(String.format("\n"));
      }
      return sb.toString();
  }
}