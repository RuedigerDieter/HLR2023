#include <stdio.h>

// Definieren Sie ein 3x3-Array Namens map, das Werte vom Typ double enthält
double map[3][3];

// Die Funktion set_temperature soll an Position [x, y] den Wert dir in das Array map eintragen
// Überprüfen Sie x und y, um mögliche Arrayüberläufe zu verhindern
void set_temperature (int x, int y, double temperature)
{
	if(x < 0 || x > 2)
	{
		printf("Fehler: x=%d liegt außerhalb des Arrays\n",x);
		return;
	}
	if(y < 0 || y > 2)
	{
		printf("Fehler: y=%d liegt außerhalb des Arrays\n",y);
		return;
	}
	map[x][y] = temperature;
}

// Die Funktion show_map soll das Array in Form einer 3x3-Matrix ausgeben
void show_map (void)
{
	for (int i = 0; i < 3; i++)
	{
		/*
		 * %.2f um die Float mit zwei Nachkommastellen darzustellen, \t um die Tabelle zu formatieren.
		 */
		for (int k = 0; k < 3; k++)
		{
			printf("%.2f \t", map[i][k]);
		}
		printf("\n");
	}
}

// Die Funktion average_value soll an Position [x, y] den Durchschnitt der 8 umgebenen
// Temperaturen in das Array map eintragen. 
// Für Werte außerhalb des Arrays soll der Wert 0 angenommen werden.
// Verwenden Sie hierfür auch die Funktion set_temperature.

/*
 * 3*3 Zellen werden addiert, für bessere Verständlichkeit wird die "Urzelle" später wieder von der Summe abgezogen, 
 * sodass der Durchschnitt ausgerechnet und auf auf die ausgewählte Zelle geschrieben werden kann.
 */
void set_to_average (int x, int y)
{
  double avg = 0;
  for (int i = x - 1; i <= x+1; i++)
  {
	if (i < 0 || i > 2)
	{
		continue;
	}
	else
	{
		for (int k = y - 1; k <= y+1; k++)
		{
			if(k < 0 || k > 2)
			{
				continue;
			}
			else
			{
				avg += map[i][k];
			}
		}
  	}
	}
	avg -= map[x][y];
	avg = avg / 8;
	set_temperature(x, y, avg);
}

// In dieser Funktion darf nichts verändert werden!
int main (void)
{
	set_temperature(0, 1, 40);
	set_temperature(1, 0, 160);
	set_temperature(1, 4, 75);
	set_temperature(1, 2, 80);
	set_temperature(2, 1, 120);

	show_map();

	set_temperature(0, 0, 20.5);
	set_temperature(0, 2, 14.8);
	set_temperature(0, 2, 22.7);
	set_temperature(2, 0, 100.2);
	set_temperature(2, 2, 20.6);
	set_temperature(2, 2, 200.2);
	set_temperature(1, 3, 201.06);
	set_temperature(1, 1, 50.5);

	show_map();
  
  set_to_average(0,0);
  set_to_average(2,0);
  set_to_average(1,2);
  
  show_map();

	return 0;
}
