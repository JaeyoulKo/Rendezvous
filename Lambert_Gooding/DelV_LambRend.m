function delV = DelV_LambRend(Mu, r1vec, v1vec, r2vec, v2vec, delT)

[~, delV1Vec, delV2Vec] = DelVVec_LambRend(Mu, r1vec, v1vec, r2vec, v2vec, delT);
delV = min(vecnorm(delV1Vec)+vecnorm(delV2Vec));

end

function [nSol, delV1Vec, delV2Vec] = DelVVec_LambRend(Mu, r1vec, v1vec, r2vec, v2vec, delT)

[nSol, v1Lam, v2Lam] = Lambert(Mu, r1vec, r2vec, delT);

delV1Vec = v1Lam-v1vec;
delV2Vec = -v2Lam+v2vec;

end

