public class AtomPair{
  Atom a1;
  Atom a2;
  public AtomPair(Atom a1, Atom a2){
    this.a1 = a1;
    this.a2 = a2;
  }
  
  public void setPair(Atom a1, Atom a2){
    this.a1 = a1;
    this.a2 = a2;
  }
  
  public Atom getA1(){
    return a1;
  }
  
  public Atom getA2(){
    return a2;
  }
  
  public String toString(){
    return String.format("[%s, %s]", a1.toString(), a2.toString());
  }
}