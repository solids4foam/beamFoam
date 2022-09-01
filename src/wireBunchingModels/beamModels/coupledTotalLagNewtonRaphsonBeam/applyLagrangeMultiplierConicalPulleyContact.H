{
    // Find active contact beam segments
    
    start = 0;
    eqn.activeConicalPulleyContacts().clear();
    
    for (label bI=0; bI<nBeams; bI++)
    {
        label nSeg = contact().splines()[bI].nSegments();

        const conicalPulleyContactListList& curPulleyContacts =
            contact().conicalPulleyContacts()[bI];

        for (label segI=0; segI<nSeg; segI++)
        {
            label globalSegI = start + segI;

            for (label pulleyI=0; pulleyI<nPulleys; pulleyI++)
            {
                const scalarField& curNormalGap =
                    curPulleyContacts[segI][pulleyI]
                   .normalGap();

                const scalarField curNormalForce =
                    mag
                    (
                        curPulleyContacts[segI][pulleyI]
                       .normalContactForce()
                    );

                for (label cI=0; cI<2; cI++)
                {
                    if (iOuterCorr()>0)
                    {
                        if
                        (
                            curNormalGap[cI] < -SMALL
                         || curNormalForce[cI] > SMALL
                        )
                        {
                            FixedList<label,3> curContact;
                            curContact[0] = bI;
                            curContact[1] = segI;
                            curContact[2] = pulleyI;

                            eqn.activeConicalPulleyContacts().insert
                            (
                                globalSegI,
                                curContact
                            );
                            break;
                        }
                    }
                }
            }
        }

        start += nSeg;
    }
}
