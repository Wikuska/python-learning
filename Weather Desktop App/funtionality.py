import requests
import os

API_KEY = os.getenv("OW_API_KEY")

def get_weather(city):
    get_weather_url = "https://api.openweathermap.org/data/2.5/weather"
    params = {
        "q": city,
        "appid": API_KEY,
    }
    response = requests.get(url = get_weather_url, params = params)
    response.raise_for_status()
    data = response.json()
    weather_id = data["weather"][0]["id"]
    weather = data["weather"][0]["description"]
    temp_in_c = round(data["main"]["temp"] - 273.15)
    return weather_id, weather, temp_in_c

