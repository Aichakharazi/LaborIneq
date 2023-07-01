function WW = solveomegas(x,aalpha,ssp, ssigmaparam, xiparam,zz, ttheta_hf, ttheta_hn, ttheta_lf, ttheta_ln, pphi_hf,  pphi_hn, pphi_lf, pphi_ln, varrho, upsilon, varepsilon, zetaparam, LY, LM, LM_hf, LM_hn, LM_lf, LM_ln, LY_hf, LY_hn, LY_lf, LY_ln, margcost)
WW = [










-x(1) + margcost * aalpha *zz * ( LY_hf* pphi_hf)^(1/varrho-1) *  (( LY_hf* pphi_hf)^(1/varrho) + (LY_hn * pphi_hn )^(1/varrho))^(varrho/ssigmaparam-1) ...
* ((( LY_hf* pphi_hf)^(1/varrho) + (LY_hn * pphi_hn )^(1/varrho))^(varrho/ssigmaparam) + ((LY_lf * pphi_lf)^(1/upsilon) + (LY_ln * pphi_ln)^(1/upsilon))^(upsilon/ssigmaparam))^(ssigmaparam-1) ...
*(LY)^((ssp-1)/ssp-1) * ( aalpha * (LY)^((ssp-1)/ssp) + (1-aalpha)* (LM)^((ssp-1)/ssp) )^(ssp/(ssp-1)-1);

-x(2) + margcost * aalpha *zz * (LY_hn * pphi_hn )^(1/varrho-1) * (( LY_hf* pphi_hf)^(1/varrho) + (LY_hn * pphi_hn )^(1/varrho))^(varrho/ssigmaparam-1) ...
* ((( LY_hf* pphi_hf)^(1/varrho) + (LY_hn * pphi_hn )^(1/varrho))^(varrho/ssigmaparam) + ((LY_lf * pphi_lf)^(1/upsilon) + (LY_ln * pphi_ln)^(1/upsilon))^(upsilon/ssigmaparam))^(ssigmaparam-1)...
*(LY)^((ssp-1)/ssp-1) * ( aalpha * (LY)^((ssp-1)/ssp) + (1-aalpha)* (LM)^((ssp-1)/ssp) )^(ssp/(ssp-1)-1);


-x(3) + margcost * aalpha *zz * (LY_lf * pphi_lf)^(1/upsilon-1) * ((LY_lf * pphi_lf)^(1/upsilon) + (LY_ln * pphi_ln)^(1/upsilon))^(upsilon/ssigmaparam-1)...
* ((( LY_hf* pphi_hf)^(1/varrho) + (LY_hn * pphi_hn )^(1/varrho))^(varrho/ssigmaparam) + ((LY_lf * pphi_lf)^(1/upsilon) + (LY_ln * pphi_ln)^(1/upsilon))^(upsilon/ssigmaparam))^(ssigmaparam-1)...
*(LY)^((ssp-1)/ssp-1) * ( aalpha * (LY)^((ssp-1)/ssp) + (1-aalpha)* (LM)^((ssp-1)/ssp) )^(ssp/(ssp-1)-1);


-x(4) + margcost * aalpha *zz * (LY_ln * pphi_ln)^(1/upsilon-1) * ((LY_lf * pphi_lf)^(1/upsilon) + (LY_ln * pphi_ln)^(1/upsilon))^(upsilon/ssigmaparam-1)...
* ((( LY_hf* pphi_hf)^(1/varrho) + (LY_hn * pphi_hn )^(1/varrho))^(varrho/ssigmaparam) + ((LY_lf * pphi_lf)^(1/upsilon) + (LY_ln * pphi_ln)^(1/upsilon))^(upsilon/ssigmaparam))^(ssigmaparam-1)...
*(LY)^((ssp-1)/ssp-1) * ( aalpha * (LY)^((ssp-1)/ssp) + (1-aalpha)* (LM)^((ssp-1)/ssp) )^(ssp/(ssp-1)-1);



-x(5) + margcost * (1-aalpha) *zz * ( LM_hf * ttheta_hf)^(1/varepsilon-1) * (( LM_hf * ttheta_hf)^(1/varepsilon) + (LM_hn * ttheta_hn)^(1/varepsilon))^(varepsilon/xiparam-1)...
* ((( LM_hf * ttheta_hf)^(1/varepsilon) + (LM_hn * ttheta_hn)^(1/varepsilon))^(varepsilon/xiparam) + ((LM_lf * ttheta_lf)^(1/zetaparam) + (LM_ln * ttheta_ln)^(1/zetaparam))^(zetaparam/xiparam))^(xiparam-1)...
*(LM)^((ssp-1)/ssp-1) * ( aalpha * (LY)^((ssp-1)/ssp) + (1-aalpha)* (LM)^((ssp-1)/ssp) )^(ssp/(ssp-1)-1);

-x(6) +  margcost * (1-aalpha) *zz * (LM_hn * ttheta_hn)^(1/varepsilon-1) * (( LM_hf * ttheta_hf)^(1/varepsilon) + (LM_hn * ttheta_hn)^(1/varepsilon))^(varepsilon/xiparam-1)...
* ((( LM_hf * ttheta_hf)^(1/varepsilon) + (LM_hn * ttheta_hn)^(1/varepsilon))^(varepsilon/xiparam) + ((LM_lf * ttheta_lf)^(1/zetaparam) + (LM_ln * ttheta_ln)^(1/zetaparam))^(zetaparam/xiparam))^(xiparam-1)...
*(LM)^((ssp-1)/ssp-1) * ( aalpha * (LY)^((ssp-1)/ssp) + (1-aalpha)* (LM)^((ssp-1)/ssp) )^(ssp/(ssp-1)-1);


-x(7) +  margcost * (1-aalpha) *zz * (LM_lf * ttheta_lf)^(1/zetaparam-1) * ((LM_lf * ttheta_lf)^(1/zetaparam) + (LM_ln * ttheta_ln)^(1/zetaparam))^(zetaparam/xiparam-1)...
* ((( LM_hf * ttheta_hf)^(1/varepsilon) + (LM_hn * ttheta_hn)^(1/varepsilon))^(varepsilon/xiparam) + ((LM_lf * ttheta_lf)^(1/zetaparam) + (LM_ln * ttheta_ln)^(1/zetaparam))^(zetaparam/xiparam))^(xiparam-1)...
*(LM)^((ssp-1)/ssp-1) * ( aalpha * (LY)^((ssp-1)/ssp) + (1-aalpha)* (LM)^((ssp-1)/ssp) )^(ssp/(ssp-1)-1);

-x(8) +  margcost * (1-aalpha) *zz * (LM_ln * ttheta_ln)^(1/zetaparam-1)  * ((LM_lf * ttheta_lf)^(1/zetaparam) + (LM_ln * ttheta_ln)^(1/zetaparam))^(zetaparam/xiparam-1)...
* ((( LM_hf * ttheta_hf)^(1/varepsilon) + (LM_hn * ttheta_hn)^(1/varepsilon))^(varepsilon/xiparam) + ((LM_lf * ttheta_lf)^(1/zetaparam) + (LM_ln * ttheta_ln)^(1/zetaparam))^(zetaparam/xiparam))^(xiparam-1)...
*(LM)^((ssp-1)/ssp-1) * ( aalpha * (LY)^((ssp-1)/ssp) + (1-aalpha)* (LM)^((ssp-1)/ssp) )^(ssp/(ssp-1)-1);
];