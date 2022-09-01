if (beam.conicalPulleys().size())
{
    Info << "Add explicit contribution of conical pulley contac force"
         << endl;

    label nPulleys = beam.conicalPulleys().size();

    label start = 0;
    label index = 0;
    
    for (label bI=0; bI<nBeams; bI++)
    {
        const conicalPulleyContactListList& curPulleyContacts =
	    beam.contact().conicalPulleyContacts()[bI];

        vectorField curq
        (
            beam.contact().splines()[bI].nSegments(),
            vector::zero
        );

	for
	(
	    label segI=0;
	    segI<beam.contact().splines()[bI].nSegments();
	    segI++
	)
	{
            for (label pulleyI=0; pulleyI<nPulleys; pulleyI++)
            {
                for (label cI=0; cI<2; cI++)
                {
                    curq[segI] +=
                        curPulleyContacts[segI][pulleyI]
                       .normalContactForce()[cI];
                }
	    }
	}

        forAll(curq, segI)
        {
            // W equation
            b(index++) -= curq[segI].x()*L[start+segI];
            b(index++) -= curq[segI].y()*L[start+segI];
            b(index++) -= curq[segI].z()*L[start+segI];
            // Theta equation
            // b(index++) -= curm[segI].x()*L[start+segI];
            // b(index++) -= curm[segI].y()*L[start+segI];
            // b(index++) -= curm[segI].z()*L[start+segI];
            index++;
            index++;
            index++;
        }

        start += curq.size();
    }

    Info << "Add implicit contribution of conical pulley contac force"
         << endl;

    start = 0;
    for (label bI=0; bI<nBeams; bI++)
    {
        const conicalPulleyContactListList& curPulleyContacts =
	    beam.contact().conicalPulleyContacts()[bI];

        label nSeg = beam.contact().splines()[bI].nSegments();

        for (label segI=0; segI<nSeg; segI++)
        {
            for (label pulleyI=0; pulleyI<nPulleys; pulleyI++)
            {
                // if
                // (
                //     pulleyI
                //  == curPulleyContacts[segI][pulleyI]
                // )
                {
                    // const scalarField& curContactDistance =
                    //     curPulleyContacts[segI][pulleyI].delta();
			
                    // const vectorField& curContactForce =
                    //     curPulleyContacts[segI][pulleyI].normalContactForce();

                    const tensorField& curContactForceDerivative =
                        curPulleyContacts[segI][pulleyI]
                       .normalContactForceDerivative();

                    // const label nSegments =
                    //     beam.contact().splines()[nbI].nSegments();

                    for (label cI=0; cI<2; cI++)
                    {
                        // Diagonal
                        label globalSegI = start + segI;
                        label ID = 6*globalSegI;

                        //- Contact force derivative (W)
                        A.coeffRef(ID,ID) +=
                            curContactForceDerivative[cI].xx()*L[globalSegI];
                        A.coeffRef(ID,ID+1) +=
                            curContactForceDerivative[cI].xy()*L[globalSegI];
                        A.coeffRef(ID,ID+2) +=
                            curContactForceDerivative[cI].xz()*L[globalSegI];

                        A.coeffRef(ID+1,ID) +=
                            curContactForceDerivative[cI].yx()*L[globalSegI];
                        A.coeffRef(ID+1,ID+1) +=
                            curContactForceDerivative[cI].yy()*L[globalSegI];
                        A.coeffRef(ID+1,ID+2) +=
                            curContactForceDerivative[cI].yz()*L[globalSegI];

                        A.coeffRef(ID+2,ID) +=
                            curContactForceDerivative[cI].zx()*L[globalSegI];
                        A.coeffRef(ID+2,ID+1) +=
                            curContactForceDerivative[cI].zy()*L[globalSegI];
                        A.coeffRef(ID+2,ID+2) +=
                            curContactForceDerivative[cI].zz()*L[globalSegI];
                    }
                }
            }
        }
            
        start += nSeg;
    }
}
