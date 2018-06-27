  Eigen::VectorXd reproject (Eigen::VectorXd streamfunction ,
                 Eigen::VectorXd eigenvector ,
                 Eigen::VectorXd transpose )  {

return  streamfunction - (transpose.cwiseProduct( streamfunction).cwiseProduct( eigenvector)); 

}   
  
    
  
  
UT_Vector3F g(GU_Detail *gdp, Eigen::VectorXd &streamfunction, GA_Offset primnum, GEO_PolyInterface &t,GEO_HedgeInterface &t2 ) {


	GEO_Hedge Hedge1 = t.polyHedge(primnum);
	GA_Offset vtx0 = t.srcPoint(Hedge1);
	GEO_Hedge  hnext = t.nextPrimitiveHedge(Hedge1);
	GA_Offset vtx1 = t.srcPoint(hnext);
	GEO_Hedge  hprev = t.prevPrimitiveHedge(Hedge1);
	GA_Offset vtx2 = t.srcPoint(hprev);
	UT_Vector3F prev = t2.hedgeVector(hprev);
	UT_Vector3F next = t2.hedgeVector(hnext);
	UT_Vector3F curr = t2.hedgeVector(Hedge1);

	UT_Vector3F	   n = cross(next, prev); n.normalize();
	fpreal area = gdp->getGEOPrimitive(primnum)->calcArea();

	UT_Vector3F  g =
		cross(n, next) *  streamfunction[vtx0] +
		cross(n, prev) *  streamfunction[vtx1] +
		cross(n, curr) *  streamfunction[vtx2];
  
g = g / (2.0 * area);
g =  cross( n, g);

return g; 
}



 Eigen::SparseMatrix<double> fvf(GU_Detail *gdp ,  Eigen::VectorXd streamfunction,Eigen::SparseMatrix<double> A_Inv) {
  int  primnums= gdp->getNumPrimitives(); 
  Eigen::MatrixXd gradreturn( primnums, 3); 
  Eigen::SparseMatrix<double> fa(gdp->getNumPoints(), gdp->getNumPoints());  
  GEO_PolyInterface t(gdp);
  GEO_HedgeInterface t2(gdp);
  GA_Offset  ptoff;
  fpreal32 sum = 0; 

      GA_FOR_ALL_PTOFF(gdp, ptoff)
      {
          GEO_Hedge h =  t.firstIncidentHedge(ptoff);

            do{ 
                        //or gdp->buildRingZeroPoints()
                        if ( t.dstPoint(h) != ptoff) {
                          
                            GEO_Hedge heq = t.nextEquivalentHedge( h);            
                            GEO_Hedge prevPrim = t.prevPrimitiveHedge(h); 
                            GEO_Hedge nextEqPrim = t.nextPrimitiveHedge(heq); 
                              // std::cout << t.dstPoint(h)<<"ta"<<ptoff <<std::endl; 
                              
                                  int hx  =t.dstPoint(h) ; 
                                  int primnum0 = t.hedgePoly(h);
                                  int primnum1 = t.hedgePoly(heq);  
                                  GEO_Primitive *prim1 = gdp->getGEOPrimitiveByIndex(primnum0);
                                  GEO_Primitive *prim2 = gdp->getGEOPrimitiveByIndex(primnum1);

                                  int src = t.srcPoint(h); 
                                  int n1 =t.srcPoint(prevPrim); 
                                  int n2 = t.dstPoint( nextEqPrim); 

                                  UT_Vector3F deriativea = gdp->getPos3( n1) - gdp->getPos3( hx);
                                  UT_Vector3F deriativeb = gdp->getPos3( hx) - gdp->getPos3( n2);
                                  UT_Vector3F vn1 = prim1->computeNormal(); 
                                  UT_Vector3F vn2 = prim2->computeNormal();    
                                  deriativea = cross( vn1, deriativea);
                                  deriativeb=  cross( vn2, deriativeb);	

                          UT_Vector3F g1 = g( gdp, streamfunction, primnum0, t, t2);
                          UT_Vector3F g2= g( gdp, streamfunction, primnum1, t, t2);
                
                              
                            fpreal32 fvd = 1.0/6.0 *  (dot( deriativea, g1) + dot(deriativeb, g2)); 
                            sum += fvd; 
                            
                            fa.coeffRef( ptoff , hx) = fvd; 
                            fa.coeffRef( ptoff , ptoff) =sum; 
                                
  

                        }
                h = t.nextIncidentHedge(h,ptoff);

            }   while (h != t.firstIncidentHedge(ptoff));
  }
  fa = A_Inv * fa; 

  return fa; 
  }
