
#let betaw = text(fill: purple)[$beta_w_z$]
#let excr = text(fill: purple)[$R_(Y ~ W | Z, U)$] 
#let excr2 = text(fill: purple)[$R^2_(Y ~ W | Z, U)$] 
#let exogr = text(fill: purple)[$R_(W ~ U)$] 
#let exogr2 = text(fill: purple)[$R^2_(W ~ U)$] 
#let betau = text(fill: blue)[$beta_u$]
#let phiu = text(fill: purple)[$phi_u$]

// = Bounds on Robustness Values for IV When Only One Identifying Assumption is Violated
//
// In this note, I derive a bounds on the robustness values for IV that are valid when the exclusion restriction is known to hold, and the only violations in the identifying assumptions come in the form of exogeneity violations. In this case, the IV robustness value will be no more than the R-squared value for the reduced-form regression of $Y$ on $W$: $"RV"_"IV" lt.eq "Cov"(Y,W)^2/("Var"(W) "Var"(Y))$
//
// - The first bound is applicable when exogeneity is known to hold, and only the exclusion restriction can be violated. In this case, the value of the IV robustness value can be no more than 
// $tau_"IV"^2 pi^(2) frac("Var"(W^(perp Z)), "Var"(Y^(perp Z)))   frac(1, 1 + 19 times R^2_(W ~ Z)) gt.eq  excr2 $ , where $pi$ is the first-stage coefficient from the first-stage regression. 


I start with the formula given in Theorem 3.1 of the original paper, which I replicate after removing references to $W^Y$: 


$
"bias"(tau_("IV")) =  underbrace(frac(1, R_(Z ~ W)) frac("sd"(W), "sd"(Z)), "Observable Scaling Factor"){
  underbrace(betaw, "Exclusion Restriction")  + betau underbrace(phiu, "Exogeneity") frac("Var"(U), "Var"(W))
}
$

I then substitute in $frac("Cov"(Z, W), "sd"(Z) "sd"(W))$ for $R_(Z~W)$ and simplify: 


$
"bias"(tau_("IV")) &=  frac("sd"(Z) "sd"(W), "Cov"(Z,W)) frac("sd"(W), "sd"(Z)) {
  betaw  + betau phiu frac("Var"(U), "Var"(W))
} \ 
"bias"(tau_("IV")) &=  frac("Var"(W), "Cov"(Z,W))  {betaw  + betau phiu frac("Var"(U), "Var"(W)) }
$

Generally, we are interested in computing the minimal violations in the IV / RV assumptions that are enough to drive the true IV effect to zero (a bias of size $tau_"IV"$). This corresponds to finding the minimum values of $exogr$, $excr$ that satisfy the following equality: 

$
abs(underbrace(frac("Cov"(W,Y), "Cov"(W,Z)), tau_
"IV"))  &=   abs( frac("Var"(W), "Cov"(Z,W)) { betaw  + betau phiu frac("Var"(U), "Var"(W))})
$


== Assuming the Exclusion Restriction Holds: 

$
abs(tau_"IV") &=   abs( frac("Var"(W), "Cov"(Z,W)) { underbrace(cancel(betaw),=0)  + betau phiu frac("Var"(U), "Var"(W))}) & ("Assume Exclusion Restiction)"\
abs(tau_"IV") &=   abs( frac("Var"(W), "Cov"(Z,W)) frac("Var"(U), "Var"(W)) betau phiu)  & ("Commute")\
abs(tau_"IV") &=   abs( frac("Var"(U), "Cov"(Z,W)) betau phiu) & ("Cancel Var"(W))\
$

I then use $phiu = R_(W ~ U) frac("sd"(W), "sd"(U))$ to write: 

$
abs(frac("Cov"(W,Y), "Cov"(W,Z))) &= abs(
   frac("Var"(U), "Cov"(Z,W)) 
   betau
   R_(W ~ U) frac("sd"(W), "sd"(U))) & 
   ("Substitute in" phiu)
   \ 
(frac("Cov"(W,Y), "Cov"(W,Z)))^2 &= (
   frac("Var"(U), "Cov"(Z,W)) 
   betau
   R_(W ~ U) frac("sd"(W), "sd"(U))
)^2  & ("Square Both Sides")
\ 
"Cov"(W,Y)^2 &= 
   "Var"(U)^2
   betau^2
   R^2_(W ~ U) frac("Var"(W), "Var"(U))
   & ("power of a product")
\ 
"Cov"(W,Y)^2 &= 
   "Var"(U) "Var"(W)
   betau^2
   R^2_(W ~ U)  
\ 
"Cov"(W,Y)^2 / ("Var"(U) "Var"(W) betau^2) &= 
   R^2_(W ~ U)   \
\ 
$

Where the right hand side is the first component of the IV robustness value. Assuming both exclusion restriction and no true effect, we have $betau^2 "Var"(U) = "Var"(Y) R^2_(Y ~ U)$. Performing this substitution, we have: 

$
\
"Cov"(W,Y)^2 / ("Var"(Y) R^2_(Y~U) "Var"(W)) &= 
   R^2_(W ~ U)  
\
$

Finally, we can take the Infimum over the left hand side to bound the RV on the right. Sending $R^2_(Y ~ U) -> 1$, we have:

$
inf("Cov"(W,Y)^2 / ("Var"(Y) R^2_(Y~U) "Var"(W))) &= inf( 
   R^2_(W ~ U))
\
"Cov"(W,Y)^2 / ("Var"(Y) "Var"(W)) &> 
   inf(R^2_(W ~ U))
\ 
R^2_(W ~ Y) &>  inf(R^2_(W ~ U))
$

The right hand side is the reduced-form R-squared, so we have just proved that the $"IV"_{R V}$ is bounded above by the reduced-form R-squared value. Importantly
however, in the last steps, we assumed that we could send $R^2_(Y ~ U) -> 1$ without violating the treatment confounding cap that $R^2_(Z ~ U | W) < .95$. This cap will not be an issue in settings where the treatment and outcome are not heavily correlated, but if the treatment and outcome are very correlated,
then it will not be possible to send $R^2_(Y ~ U) -> 1$ without also breaking the treatment confounding cap. Specifically, this cap will never be an issue if $R^2_(Z ~ Y | W) < .95)$.



