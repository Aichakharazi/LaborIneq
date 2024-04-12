function WW = solveomegas(x, LY, LM, aalpha,ssp, ssigmaparam, xiparam,zz, theta_hf, theta_hn, theta_lf, theta_ln, phi_hf, phi_hn, phi_lf, phi_ln, varrho, upsilon, varepsilon, zetaparam, l_hftm, l_hntm, l_lftm, l_lntm, l_hfty, l_hnty, l_lfty, l_lnty, margcost)
WW = [
%  LY, LM,
% l_hftm, l_hntm, l_lftm, l_lntm, l_hfty, l_hnty, l_lfty, l_lnty

 



-x(1) + margcost * aalpha *zz* phi_hf * ( l_hfty)^(1/varrho-1) *  ( (l_hfty)^(1/varrho) * phi_hf + (l_hnty)^(1/varrho) * phi_hn )^(varrho/ssigmaparam-1) ...
* ((( l_hfty)^(1/varrho) * phi_hf + (l_hnty)^(1/varrho) * phi_hn )^(varrho/ssigmaparam) + ((l_lfty)^(1/upsilon)  * phi_lf+ (l_lnty)^(1/upsilon) * phi_ln)^(upsilon/ssigmaparam))^(ssigmaparam-1) ...
*(LY)^(1/ssp-1) * ( aalpha * (LY)^(1/ssp) + (1-aalpha)* (LM)^(1/ssp) )^(ssp-1);



-x(2) + margcost * aalpha *zz * phi_hn  * (l_hnty )^(1/varrho-1) * (( l_hfty)^(1/varrho)* phi_hf + (l_hnty)^(1/varrho) * phi_hn )^(varrho/ssigmaparam-1) ...
* ((( l_hfty)^(1/varrho) * phi_hf + (l_hnty)^(1/varrho) * phi_hn )^(varrho/ssigmaparam) + ((l_lfty)^(1/upsilon) * phi_lf + (l_lnty)^(1/upsilon) * phi_ln)^(upsilon/ssigmaparam))^(ssigmaparam-1)...
*(LY)^(1/ssp-1) * ( aalpha * (LY)^(1/ssp) + (1-aalpha)* (LM)^(1/ssp) )^(ssp-1);


-x(3) + margcost * aalpha *zz * phi_lf* (l_lfty)^(1/upsilon-1)  *  ((l_lfty)^(1/upsilon) * phi_lf + (l_lnty)^(1/upsilon) * phi_ln)^(upsilon/ssigmaparam-1)...
* ((( l_hfty)^(1/varrho) * phi_hf + (l_hnty)^(1/varrho) * phi_hn )^(varrho/ssigmaparam) + ((l_lfty)^(1/upsilon)  * phi_lf+ (l_lnty)^(1/upsilon) * phi_ln)^(upsilon/ssigmaparam))^(ssigmaparam-1)...
*(LY)^(1/ssp-1) * ( aalpha * (LY)^(1/ssp) + (1-aalpha)* (LM)^(1/ssp) )^(ssp-1);


-x(4) + margcost * aalpha *zz  * phi_ln * (l_lnty )^(1/upsilon-1) * ((l_lfty)^(1/upsilon) * phi_lf + (l_lnty)^(1/upsilon) * phi_ln)^(upsilon/ssigmaparam-1)...
* ((( l_hfty)^(1/varrho) * phi_hf + (l_hnty )^(1/varrho) * phi_hn)^(varrho/ssigmaparam) + ((l_lfty)^(1/upsilon) * phi_lf + (l_lnty)^(1/upsilon) * phi_ln)^(upsilon/ssigmaparam))^(ssigmaparam-1)...
*(LY)^(1/ssp-1) * ( aalpha * (LY)^(1/ssp) + (1-aalpha)* (LM)^(1/ssp) )^(ssp-1);



-x(5) + margcost * (1-aalpha) *zz * theta_hf * ( l_hftm )^(1/varepsilon-1) * (( l_hftm)^(1/varepsilon) * theta_hf + (l_hntm)^(1/varepsilon) * theta_hn)^(varepsilon/xiparam-1)...
* ((( l_hftm)^(1/varepsilon)  * theta_hf+ (l_hntm)^(1/varepsilon) * theta_hn)^(varepsilon/xiparam) + ((l_lftm)^(1/zetaparam) * theta_lf + (l_lntm)^(1/zetaparam) * theta_ln)^(zetaparam/xiparam))^(xiparam-1)...
*(LM)^(1/ssp-1) * ( aalpha * (LY)^(1/ssp) + (1-aalpha)* (LM)^(1/ssp) )^(ssp-1);

-x(6) +  margcost * (1-aalpha) *zz * theta_hn * (l_hntm )^(1/varepsilon-1) * (( l_hftm)^(1/varepsilon)  * theta_hf+ (l_hntm)^(1/varepsilon) * theta_hn)^(varepsilon/xiparam-1)...
* ((( l_hftm)^(1/varepsilon) * theta_hf + (l_hntm)^(1/varepsilon) * theta_hn)^(varepsilon/xiparam) + ((l_lftm)^(1/zetaparam) * theta_lf + (l_lntm)^(1/zetaparam) * theta_ln)^(zetaparam/xiparam))^(xiparam-1)...
*(LM)^(1/ssp-1) * ( aalpha * (LY)^(1/ssp) + (1-aalpha)* (LM)^(1/ssp) )^(ssp-1);


-x(7) +  margcost * (1-aalpha) *zz * theta_lf * (l_lftm )^(1/zetaparam-1) * ((l_lftm)^(1/zetaparam) * theta_lf + (l_lntm)^(1/zetaparam) * theta_ln)^(zetaparam/xiparam-1)...
* ((( l_hftm)^(1/varepsilon)  * theta_hf+ (l_hntm)^(1/varepsilon) * theta_hn)^(varepsilon/xiparam) + ((l_lftm)^(1/zetaparam)  * theta_lf + (l_lntm)^(1/zetaparam) * theta_ln)^(zetaparam/xiparam))^(xiparam-1)...
*(LM)^(1/ssp-1) * ( aalpha * (LY)^(1/ssp) + (1-aalpha)* (LM)^(1/ssp) )^(ssp-1);

-x(8) +  margcost * (1-aalpha) *zz * theta_ln * (l_lntm )^(1/zetaparam-1)  * ((l_lftm)^(1/zetaparam) * theta_lf + (l_lntm)^(1/zetaparam) * theta_ln)^(zetaparam/xiparam-1)...
* ((( l_hftm)^(1/varepsilon)  * theta_hf+ (l_hntm)^(1/varepsilon) * theta_hn)^(varepsilon/xiparam) + ((l_lftm)^(1/zetaparam)  * theta_lf + (l_lntm)^(1/zetaparam) * theta_ln )^(zetaparam/xiparam))^(xiparam-1)...
*(LM)^(1/ssp-1) * ( aalpha * (LY)^(1/ssp) + (1-aalpha)* (LM)^(1/ssp) )^(ssp-1);


];