function v = OP_DIAGSPLIT()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = casadiMEX(0, 82);
  end
  v = vInitialized;
end
