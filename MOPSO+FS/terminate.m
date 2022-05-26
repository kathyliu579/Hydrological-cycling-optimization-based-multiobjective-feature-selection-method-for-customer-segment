function y=terminate(evaluation,maxEvaluation)
  if evaluation>maxEvaluation
       y=false;
       return;
  end
  y=true;
end