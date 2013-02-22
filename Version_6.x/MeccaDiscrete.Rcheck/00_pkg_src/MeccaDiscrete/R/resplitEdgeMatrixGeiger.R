resplitEdgeMatrixGeiger <-
function (z, phy, node, cutAtStem=T, n2=1) 
{
	newTag<-max(z[,7])+1
    rootnode <- length(phy$tip.label) + 1

	dec2<-numeric(length(z[,2]))
	for(i in 1:max(z[,2])) {
		ss<-z[,2]==i
		nn<-sum(ss)
		oo<-order(z[ss,4], decreasing=T)
		dec2[ss]<-oo
		}


    if (node >= rootnode) {
        node.desc <- node
        pos <- 1
        
        if(cutAtStem) {
        	row1<-which(z[, 2] == node.desc[1] & dec2 >= n2)
        	row2<-which(z[, 1] == node.desc[1] & dec2 >= n2)
 			row<-c(row1, row2)
 		} else {
 			row<-which(z[, 2] == node.desc[1] & dec2 >= n2)
 			}       
        base<-min(z[row,7])
        ok<- z[row,7]==base 
        if(sum(ok)>0)
               z[row[ok],7] <- newTag
           

        while (pos != (length(node.desc) + 1)) {
            temp <- node.sons(phy, node.desc[pos])
            temp <- temp[temp > rootnode]
            if(length(temp)!=0) {
             for (k in 1:length(temp)) {
             	row<-which(z[, 1] == temp[k])
                ok<- z[row,7]==base 
                if(sum(ok)>0)
                	z[row[ok],7] <- newTag
            	}
             node.desc <- c(node.desc, temp)
            }
            pos <- pos + 1
        }
    } else if (node > 0) {
        ok<-z[,2]==node & dec2>=n2
        z[ok,7] <- newTag
    }
    return(z)
}
