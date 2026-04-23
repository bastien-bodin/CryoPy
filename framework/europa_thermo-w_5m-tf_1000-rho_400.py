width = 5.0
tf_days = 1000.0
dp = 0.05
rho_snow = 400.0

path2dir = '../outputs/' # You can also set as absolute path
prefix_output = f'EuropaThermoRad'
saving_directory = f'{path2dir}{prefix_output}-w_{int(width)}-rho_{rho_snow}-tf_{tf_days}_output/'

def main():
    from apps.LayerCool import ThermoSnowEuropa

    MyApp = ThermoSnowEuropa(
        fname=prefix_output, output_dir=saving_directory,
        dp=dp,rho_snow=rho_snow,W=width, tf_days=tf_days
    )
    MyApp.run()
    return None

if __name__=="__main__":
    main()
