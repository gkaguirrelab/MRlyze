function writeme(CO,n,fname,HRF,LRF)
Hstring = [ ';VB98'; ';REF1'; ';    '; ];
ORs = sprintf('; Cardiac Rate: %f',HRF);
Ms = sprintf('; File length: %u',length(CO));
writetotext(CO(:,1),[sprintf('H-Sin1-(%u)-',n) fname],Hstring,ORs,Ms)
writetotext(CO(:,2),[sprintf('H-Cos1-(%u)-',n) fname],Hstring,ORs,Ms)
writetotext(CO(:,3),[sprintf('H-Sin2-(%u)-',n) fname],Hstring,ORs,Ms)
writetotext(CO(:,4),[sprintf('H-Cos2-(%u)-',n) fname],Hstring,ORs,Ms)
ORs = sprintf('; Resp Rate: %f',LRF);
writetotext(CO(:,5),[sprintf('L-Sin1-(%u)-',n) fname],Hstring,ORs,Ms)
writetotext(CO(:,6),[sprintf('L-Cos1-(%u)-',n) fname],Hstring,ORs,Ms)
writetotext(CO(:,7),[sprintf('L-Sin2-(%u)-',n) fname],Hstring,ORs,Ms)
writetotext(CO(:,8),[sprintf('L-Cos2-(%u)-',n) fname],Hstring,ORs,Ms)