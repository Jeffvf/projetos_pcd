public class TrafficController {
    public Boolean getLeft, getRight;
    public long left, right;

    public TrafficController(){
        left = 0;
        right = 0;
        getRight = false;
        getLeft = false;
    }
    
    public void enterLeft(long ID){
        try{
        	Thread.sleep(1);
        	} catch (Exception e){}
        
        while(true){
            getLeft = true;
            left = ID;
            
            if(right != 0 && right != ID){
                getLeft = false;
                
                continue;
            }
            
            right = ID;
            
            if(left != ID){
                getLeft = false;
                
                while(getRight);{
	                if(right != ID) 
	                	continue;
	                else{
	                    getLeft = true;
	                    
	                    break;
	                }
                }
            }else 
            	break;
        
        }
    }
    
    public void enterRight(long ID){
        while(true){
            try{
            	Thread.sleep(1);
            	} catch (Exception e){}
            
            left = ID;
            getRight = true;
            
            if(right != 0){
                getRight = false;
                
                continue;
            }
            
            right = ID;
            
            if(left != ID){
                getRight = false;
                
                while(getLeft);{
	                if(right != ID) 
	                	continue;
	                else{
	                    getRight = true;
	                    
	                    break;
	                }
                }
            }else 
            	break;
        }
    }
    
    public void leaveLeft(){
        right = 0;
        
        getLeft = false;
    }
    
    public void leaveRight(){
        right = 0;
        
        getRight = false;
    }

}