# Define placeholders for the variables
fuel_pin = 'FUEL'
moderator_pin = 'MODERATOR'

# Your list structure with placeholders
rings_1 = [
    [moderator_pin, fuel_pin, fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, fuel_pin] * 6,
    [fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, moderator_pin, fuel_pin] * 6,
    [moderator_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin] * 6,
    [fuel_pin, fuel_pin, moderator_pin, fuel_pin, fuel_pin, fuel_pin, moderator_pin, fuel_pin] * 6,
    [moderator_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin] * 6,
    [fuel_pin, fuel_pin, moderator_pin, fuel_pin, moderator_pin, fuel_pin] * 6,
    [moderator_pin, fuel_pin, fuel_pin, fuel_pin, fuel_pin] * 6,
    [fuel_pin, fuel_pin, moderator_pin, fuel_pin] * 6,
    [moderator_pin, fuel_pin, fuel_pin] * 6,
    [fuel_pin, moderator_pin] * 6,
    [fuel_pin] * 6,
    [moderator_pin]
]