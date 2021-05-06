template<typename T>
std::vector<double> linspace(T start_in, T end_in, int num_in)
{
  std::vector<double> linspaced;

  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);

  if (num == 0) { return linspaced; }
  if (num == 1)
    {
      linspaced.push_back(start);
      return linspaced;
    }

  double delta = (end - start) / (num - 1);

  for(int i=0; i < num-1; ++i)
    {
      linspaced.push_back(start + delta * i);
    }
  linspaced.push_back(end); // I want to ensure that start and end
                            // are exactly the same as the input
  return linspaced;
}

  double convertToRadians(double angle)
  {
    double pi = 2*acos(0.0);
    double output = (angle / 360)*2*pi;
    return output;
  }

  double applyPeriodicBC(double num, double length)
  {
    double bound = length;
    if (num < 0)
    {
      return num + bound;
    }
    else if (num > bound) {
      return num - bound;
    } else {
      return num;
    }

  }
