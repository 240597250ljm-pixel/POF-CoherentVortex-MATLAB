function s = nearFieldStats(Pi,near,background)
% Each spatial pixel counted once, even if covered by multiple vortices.
% C_r=E[ (-Pi)*I(Pi<0) | region r] = p_r*m_r .
assert(isequal(size(Pi),size(near),size(background)));
assert(~any(near(:)&background(:)),'Masks overlap.');
a=Pi(near & isfinite(Pi));b=Pi(background & isfinite(Pi));
s=struct('nearCount',numel(a),'backgroundCount',numel(b),...
    'activityRatio',NaN,'negativeContributionRatio',NaN,...
    'negativeFrequencyRatio',NaN,'negativeConditionalMagnitudeRatio',NaN,...
    'identityError',NaN);
if isempty(a)||isempty(b),return;end
s.activityRatio=safeDiv(mean(abs(a)),mean(abs(b)));
pA=mean(a<0);pB=mean(b<0);
s.negativeFrequencyRatio=safeDiv(pA,pB);
if any(a<0)&&any(b<0)
    mA=mean(-a(a<0));mB=mean(-b(b<0));
    s.negativeConditionalMagnitudeRatio=safeDiv(mA,mB);
end
cA=mean(max(-a,0));cB=mean(max(-b,0));
s.negativeContributionRatio=safeDiv(cA,cB);
if isfinite(s.negativeContributionRatio) && isfinite(s.negativeFrequencyRatio) && ...
      isfinite(s.negativeConditionalMagnitudeRatio)
    s.identityError=s.negativeContributionRatio - ...
       s.negativeFrequencyRatio*s.negativeConditionalMagnitudeRatio;
end
end
function y=safeDiv(a,b)
if abs(b)<=eps||~isfinite(b),y=NaN;else,y=a/b;end
end
